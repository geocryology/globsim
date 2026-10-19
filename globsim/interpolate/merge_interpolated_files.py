"""
Print the NCO / nccopy commands needed to merge globsim interpolated files
from disjoint time periods into one file per kind (sa, sf, pl, pl_sur), chunked
for reading one station at a time over long time spans (i.e. for scaling).

Nothing is executed: the script only reads file headers and time coordinates,
then prints the commands. The pipeline for each kind is:

  1. ncks    : make `time` a record dimension (and trim overlapping records)
  2. ncrcat  : concatenate along time (uncompressed, keeps input chunking)
  3. nccopy  : rechunk to (time=long, level=all, station=1), make `time` fixed
  4. ncatted : tag the file with chunk_strategy=scaling
  5. verify, then delete intermediates

Usage:
    python -m globsim.merge_interpolated_files DIR_1950_2011 DIR_2011_2025 -o OUTDIR
    python -m globsim.merge_interpolated_files a_sa.nc b_sa.nc -o OUTDIR --memory-gb 16
"""
import argparse
import math
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import netCDF4 as nc
import numpy as np

from globsim.chunking import MIN_CHK_BYTES, SCL_CHK_BYTES

# {reanalysis}_{kind}_{list_name}.nc  and  {reanalysis}_pl_{list_name}_surface.nc
FILE_PATTERN = re.compile(r"^(?P<prefix>.+?)_(?P<kind>sa|sf|pl|to|sc)_(?P<list_name>.+?)(?P<surface>_surface)?\.nc$")
SKIP_SUFFIXES = ("_original.nc", "_chk.nc", "_rec.nc")  # leftovers from rechunking / previous merges
INVARIANT_KINDS = ("to", "sc")
KIND_ORDER = ("sa", "sf", "pl", "pl_sur", "to", "sc")

TIME, STATION, LEVEL, NUMBER = "time", "station", "level", "number"


@dataclass
class FileInfo:
    path: Path
    nt: int
    nl: int  # 0 if no level dimension
    ns: int
    nnum: int  # 1 if no ensemble dimension
    time_is_record: bool
    time_units: str
    calendar: str
    times: np.ndarray
    chunking: Optional[dict]
    stations: np.ndarray
    levels: Optional[np.ndarray]
    data_vars: set
    bytes_per_element: int
    data_bytes: int  # uncompressed size of all (time, ..., station) variables
    interpolation_complete: Optional[bool]


# ---------------------------------------------------------------------------
# Finding files
# ---------------------------------------------------------------------------

def classify_file(path: Path) -> Optional[tuple]:
    """Return (prefix, kind) for a globsim interpolated file name, e.g. ('era5', 'pl_sur'), or None."""
    m = FILE_PATTERN.match(path.name)
    if m is None:
        return None
    kind = m["kind"]
    if m["surface"]:
        if kind != "pl":
            return None
        kind = "pl_sur"
    return m["prefix"], kind


def find_input_files(inputs: list[Path], kinds: Optional[list[str]] = None) -> dict[tuple, list[Path]]:
    """
    Group the input files that should be merged together.

    inputs: either 2+ files or 2+ directories (one per time period).
    For directories, each directory must contain one file of each (reanalysis, kind) to be merged.
    For files, files are grouped by (reanalysis, kind); unrecognized names form a single group.
    Returns {(prefix, kind): [path, ...]} with one path per time period (unordered).
    """
    if all(p.is_dir() for p in inputs):
        groups: dict[tuple, list[Path]] = {}
        for d in inputs:
            for f in sorted(d.glob("*.nc")):
                if f.name.endswith(SKIP_SUFFIXES):
                    continue
                key = classify_file(f)
                if key is None:
                    continue
                if any(p.parent == d for p in groups.get(key, [])):
                    raise ValueError(f"More than one {key[1]} file for '{key[0]}' in {d}: "
                                     f"{[p.name for p in groups[key] if p.parent == d] + [f.name]}")
                groups.setdefault(key, []).append(f)

        for key, files in groups.items():
            missing = [str(d) for d in inputs if not any(f.parent == d for f in files)]
            if missing:
                print(f"WARNING: no {key[0]} {key[1]} file in {missing}; it will only be merged from the other directories",
                      file=sys.stderr)

    elif all(p.is_file() for p in inputs):
        groups = {}
        for f in inputs:
            key = classify_file(f) or ("merged", "unknown")
            groups.setdefault(key, []).append(f)
    else:
        raise ValueError("Provide either 2+ files or 2+ directories, not a mix (or a path that doesn't exist).")

    if kinds:
        groups = {k: v for k, v in groups.items() if k[1] in kinds or k[1] == "unknown"}

    return dict(sorted(groups.items(), key=lambda kv: (kv[0][0], KIND_ORDER.index(kv[0][1]) if kv[0][1] in KIND_ORDER else 99)))


def read_file_info(path: Path) -> FileInfo:
    """Read dimensions, time coordinate, chunking and station metadata from an interpolated file."""
    with nc.Dataset(path) as ncf:
        dims = ncf.dimensions
        if TIME not in dims or STATION not in dims:
            raise ValueError(f"{path} has no time/station dimensions")

        data_vars = {name for name, v in ncf.variables.items()
                     if TIME in v.dimensions and STATION in v.dimensions}
        chunking, bpe, data_bytes = None, 0, 0
        for name in sorted(data_vars):
            v = ncf.variables[name]
            bpe = max(bpe, v.dtype.itemsize)
            data_bytes += int(np.prod(v.shape)) * v.dtype.itemsize
            chk = v.chunking()
            if chunking is None and chk != "contiguous":
                chunking = dict(zip(v.dimensions, chk))

        tvar = ncf.variables[TIME]
        complete = getattr(ncf, "globsim_interpolate_success", None)

        return FileInfo(
            path=path,
            nt=dims[TIME].size,
            nl=dims[LEVEL].size if LEVEL in dims else 0,
            ns=dims[STATION].size,
            nnum=dims[NUMBER].size if NUMBER in dims else 1,
            time_is_record=dims[TIME].isunlimited(),
            time_units=tvar.units,
            calendar=getattr(tvar, "calendar", "standard"),
            times=np.asarray(tvar[:], dtype="f8"),
            chunking=chunking,
            stations=np.asarray(ncf.variables[STATION][:]) if STATION in ncf.variables else np.arange(dims[STATION].size),
            levels=np.asarray(ncf.variables[LEVEL][:]) if LEVEL in ncf.variables else None,
            data_vars=data_vars,
            bytes_per_element=bpe or 4,
            data_bytes=data_bytes,
            interpolation_complete=None if complete is None else bool(int(complete)),
        )


# ---------------------------------------------------------------------------
# Ordering, overlap and compatibility
# ---------------------------------------------------------------------------

def check_compatible(infos: list[FileInfo]) -> list[str]:
    """Return a list of problems that would make the files impossible to concatenate along time."""
    ref, errors = infos[0], []
    for info in infos[1:]:
        name = str(info.path)
        if info.ns != ref.ns or not np.array_equal(info.stations, ref.stations):
            errors.append(f"{name}: stations differ from {ref.path} (count or order)")
        if info.nl != ref.nl or (ref.levels is not None and not np.array_equal(info.levels, ref.levels)):
            errors.append(f"{name}: pressure levels differ from {ref.path}")
        if info.nnum != ref.nnum:
            errors.append(f"{name}: ensemble 'number' dimension differs from {ref.path}")
        if info.data_vars != ref.data_vars:
            errors.append(f"{name}: variables differ from {ref.path}: "
                          f"{sorted(info.data_vars ^ ref.data_vars)}")
        if info.calendar != ref.calendar:
            errors.append(f"{name}: calendar '{info.calendar}' differs from '{ref.calendar}'")
    return errors


def order_and_trim(infos: list[FileInfo]) -> tuple[list[tuple[FileInfo, int]], list[str]]:
    """
    Sort files by start time and work out how many leading records of each file
    must be dropped so that it starts strictly after the previous file ends.

    Returns ([(info, first_index_to_keep), ...], warnings)
    """
    ref = infos[0]

    def to_ref_units(info: FileInfo) -> np.ndarray:
        if info.time_units == ref.time_units:
            return info.times
        dates = nc.num2date(info.times, info.time_units, info.calendar)
        return np.asarray(nc.date2num(dates, ref.time_units, ref.calendar), dtype="f8")

    t = {id(i): to_ref_units(i) for i in infos}
    ordered = sorted(infos, key=lambda i: t[id(i)][0])
    warnings, result = [], []
    prev_last, step = None, None

    for info in ordered:
        times = t[id(info)]
        if len(times) > 1:
            if np.any(np.diff(times) <= 0):
                warnings.append(f"{info.path.name}: time is not strictly increasing")
            step = step or float(np.median(np.diff(times)))
        start = 0
        if prev_last is not None:
            start = int(np.searchsorted(times, prev_last, side="right"))
            if start >= len(times):
                raise ValueError(f"{info.path.name} lies entirely within the previous file's time range")
            if start > 0:
                warnings.append(f"{info.path.name}: dropping {start} overlapping leading record(s) "
                                f"(<= {fmt_time(prev_last, ref)})")
            gap = times[start] - prev_last
            if step and gap > 1.5 * step:
                warnings.append(f"{info.path.name}: gap of {gap / step - 1:.0f} missing time step(s) between "
                                f"{fmt_time(prev_last, ref)} and {fmt_time(times[start], ref)}")
        result.append((info, start))
        prev_last = times[-1]

    return result, warnings


def fmt_time(value: float, ref: FileInfo) -> str:
    return str(nc.num2date(value, ref.time_units, ref.calendar))


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------

def calculate_merged_chunking(nt: int, nl: int, nnum: int = 1, bytes_per_element: int = 4,
                              target_bytes: int = SCL_CHK_BYTES) -> dict:
    """
    Optimal chunking for reading one station at a time over long time spans.

    station = 1, level = all, number = 1, and time as long as fits in target_bytes.
    The time chunk is then balanced so the series splits into equal chunks
    (no tiny remainder chunk at the end).
    """
    c_lev = max(nl, 1)
    c_time = max(1, min(nt, target_bytes // (bytes_per_element * c_lev)))
    n_chunks = math.ceil(nt / c_time)
    c_time = math.ceil(nt / n_chunks)

    chunks = {TIME: c_time}
    if nnum > 1:
        chunks[NUMBER] = 1
    if nl:
        chunks[LEVEL] = nl
    chunks[STATION] = 1
    return chunks


def intermediate_chunking(info: FileInfo, nt_total: int) -> dict:
    """
    Chunking for the record-dimension / concatenated intermediates: keep the input
    file's chunk shape so ncks/ncrcat are close to straight copies, but avoid tiny
    record chunks (an NCO default can be 1 record per chunk, which is very slow).
    """
    existing = info.chunking or {}
    c_stat = existing.get(STATION, min(info.ns, 16))
    c_lev = existing.get(LEVEL, max(info.nl, 1))
    c_num = existing.get(NUMBER, 1)
    min_time = MIN_CHK_BYTES // (info.bytes_per_element * c_lev * c_stat * c_num)
    c_time = min(nt_total, max(existing.get(TIME, 0), min_time, 1))

    chunks = {TIME: c_time}
    if info.nnum > 1:
        chunks[NUMBER] = c_num
    if info.nl:
        chunks[LEVEL] = c_lev
    chunks[STATION] = c_stat
    return chunks


def nccopy_memory_settings(chunks: dict, nt: int, nl: int, ns: int, nnum: int,
                           bytes_per_element: int, memory_bytes: int) -> tuple[int, int, int, list[str]]:
    """
    Choose nccopy's copy buffer (-m), output chunk cache (-h) and cache slots (-e).

    Best case: the largest variable fits in -m, so each variable is read and written once.
    Otherwise: -h must hold one full time-band of output chunks (ns * nnum chunks) so that
    each output chunk stays in cache until it's complete instead of being re-compressed
    and re-written for every slab nccopy reads.
    """
    notes = []
    var_bytes = nt * max(nl, 1) * ns * nnum * bytes_per_element
    chunk_bytes = chunks[TIME] * max(nl, 1) * bytes_per_element
    band_bytes = ns * nnum * chunk_bytes
    n_chunks = ns * nnum * math.ceil(nt / chunks[TIME])

    if var_bytes * 1.1 <= memory_bytes:
        m = var_bytes
        h = min(max(64 * 1024**2, chunk_bytes * 4), memory_bytes - m) if memory_bytes > m else chunk_bytes
    else:
        h = int(band_bytes * 1.1)
        if h > memory_bytes * 0.75:
            notes.append(f"one time-band of output chunks needs {human(band_bytes)} of chunk cache but only "
                         f"{human(memory_bytes)} was allowed; chunks will be rewritten repeatedly. "
                         f"Increase --memory-gb or decrease --target-mb.")
            h = int(memory_bytes * 0.75)
        m = max(memory_bytes - h, chunk_bytes)
        notes.append(f"largest variable ({human(var_bytes)}) doesn't fit in memory; copying in slabs with a "
                     f"{human(h)} chunk cache")

    e = next_prime(max(1009, 2 * min(n_chunks, max(1, h // chunk_bytes))))
    return m, h, e, notes


def next_prime(n: int) -> int:
    def is_prime(k):
        return k > 1 and all(k % p for p in range(2, int(k ** 0.5) + 1))
    while not is_prime(n):
        n += 1
    return n


def human(n_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n_bytes < 1024 or unit == "TB":
            return f"{n_bytes:.1f} {unit}" if unit != "B" else f"{int(n_bytes)} B"
        n_bytes /= 1024


def nccopy_size(n_bytes: int) -> str:
    """nccopy accepts K/M/G suffixes; round up to whole MB."""
    mb = math.ceil(n_bytes / 1024**2)
    return f"{mb // 1024}G" if mb % 1024 == 0 else f"{mb}M"


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def nco_chunk_args(chunks: dict) -> list[str]:
    args = ["--cnk_plc=all", "--cnk_map=dmn"]
    for dim, size in chunks.items():
        args += ["--cnk_dmn", f"{dim},{size}"]
    return args


def record_dimension_command(info: FileInfo, start: int, chunks: dict, tmp_dir: Path) -> tuple[Optional[list[str]], Path]:
    """
    Command that makes `time` a record dimension (required by ncrcat) and drops
    `start` overlapping leading records. Returns (command or None, file to concatenate).
    """
    if info.time_is_record and start == 0:
        return None, info.path

    out = tmp_dir / info.path.name.replace(".nc", f"_{info.path.parent.name}_rec.nc")
    cmd = ["ncks", "-O", "-4", "-L", "0"]
    if not info.time_is_record:
        cmd += ["--mk_rec_dmn", TIME]
    if start:
        cmd += ["-d", f"{TIME},{start},"]
    cmd += nco_chunk_args(chunks) + [str(info.path), str(out)]
    return cmd, out


def merge_command(paths: list[Path], chunks: dict, out: Path) -> list[str]:
    """ncrcat along the record dimension. Intermediate is left uncompressed for speed."""
    return ["ncrcat", "-O", "-4", "-L", "0"] + nco_chunk_args(chunks) + [str(p) for p in paths] + [str(out)]


def rechunk_command(src: Path, out: Path, chunks: dict, m: int, h: int, e: int, deflate: int) -> list[str]:
    spec = ",".join(f"{dim}/{size}" for dim, size in chunks.items())
    cmd = ["nccopy", "-k", "nc4", "-u"]
    if deflate:
        cmd += ["-d", str(deflate)]
    return cmd + ["-c", spec, "-m", nccopy_size(m), "-h", nccopy_size(h), "-e", str(e), str(src), str(out)]


def format_command(cmd: list[str], shell: str) -> str:
    if shell == "powershell":
        return " ".join(a if re.fullmatch(r"[\w\-./:\\,=]+", a) else "'" + a.replace("'", "''") + "'" for a in cmd)
    return shlex.join(cmd)


def remove_command(paths: list[Path], shell: str) -> str:
    if shell == "powershell":
        return "Remove-Item " + ", ".join(format_command([str(p)], shell) for p in paths)
    return shlex.join(["rm"] + [str(p) for p in paths])


def plan_group(key: tuple, paths: list[Path], out_dir: Path, tmp_dir: Path, args) -> list[str]:
    """Build the commented list of commands to merge one group of files."""
    prefix, kind = key
    lines = [f"# {'=' * 76}", f"# {prefix} {kind}: {len(paths)} files", f"# {'=' * 76}"]

    if kind in INVARIANT_KINDS:
        lines.append(f"# invariant file (no time dimension): no merge needed, copy one of them")
        cp = ["Copy-Item"] if args.shell == "powershell" else ["cp"]
        lines.append(format_command(cp + [str(paths[0]), str(out_dir / paths[0].name)], args.shell))
        return lines

    if len(paths) < 2:
        lines.append(f"# only one file ({paths[0]}); nothing to merge. Skipping.")
        return lines

    infos = [read_file_info(p) for p in paths]
    errors = check_compatible(infos)
    if errors:
        lines += [f"# ERROR: {e}" for e in errors] + ["# Skipping this group."]
        return lines

    ordered, warnings = order_and_trim(infos)
    ref = ordered[0][0]
    nt = sum(i.nt - s for i, s in ordered)
    out = out_dir / ref.path.name
    if any(out.resolve() == p.resolve() for p in paths):
        lines.append(f"# ERROR: output {out} would overwrite an input file; choose another output directory.")
        return lines

    for info, start in ordered:
        last = fmt_time(info.times[-1], info)
        first = fmt_time(info.times[start], info)
        lines.append(f"#   {info.path}  [{first} .. {last}]  nt={info.nt - start}"
                     f"{'' if info.time_is_record else '  (time is fixed)'}")
        if info.interpolation_complete is False:
            warnings.append(f"{info.path.name}: globsim_interpolate_success is 0; interpolation may be incomplete")
    if len({i.time_units for i in infos}) > 1:
        warnings.append("time units differ between files; ncrcat rebases them to the first file's units — "
                        "check the merged time axis")
    lines += [f"# WARNING: {w}" for w in warnings]

    final = calculate_merged_chunking(nt, ref.nl, ref.nnum, ref.bytes_per_element, int(args.target_mb * 1024**2))
    inter = intermediate_chunking(ref, nt)
    m, h, e, notes = nccopy_memory_settings(final, nt, ref.nl, ref.ns, ref.nnum, ref.bytes_per_element,
                                            int(args.memory_gb * 1024**3))
    n_time_chunks = math.ceil(nt / final[TIME])
    data_bytes = sum(i.data_bytes * (i.nt - s) / i.nt for i, s in ordered)

    lines += [
        f"# merged: nt={nt}, levels={ref.nl or '-'}, stations={ref.ns}"
        f"{f', members={ref.nnum}' if ref.nnum > 1 else ''}, ~{human(data_bytes)} uncompressed",
        f"# final chunks: {final}  (~{human(final[TIME] * max(ref.nl, 1) * ref.bytes_per_element)} each, "
        f"{n_time_chunks} chunk(s) per station read)",
        f"# intermediate chunks: {inter}",
    ] + [f"# NOTE: {n}" for n in notes]

    lines.append("\n# 1. make time a record dimension (and drop overlapping records)")
    to_merge, intermediates = [], []
    for info, start in ordered:
        cmd, path = record_dimension_command(info, start, inter, tmp_dir)
        to_merge.append(path)
        if cmd is None:
            lines.append(f"#    {info.path.name}: time is already a record dimension; used as-is")
        else:
            intermediates.append(path)
            lines.append(format_command(cmd, args.shell))

    merged_tmp = tmp_dir / ref.path.name.replace(".nc", "_merged_rec.nc")
    intermediates.append(merged_tmp)
    lines.append("\n# 2. concatenate along time")
    lines.append(format_command(merge_command(to_merge, inter, merged_tmp), args.shell))

    lines.append("\n# 3. rechunk for single-station reads (and make time a fixed dimension again)")
    lines.append(format_command(rechunk_command(merged_tmp, out, final, m, h, e, args.deflate), args.shell))

    lines.append("\n# 4. tag as chunked for scaling (matches globsim.chunking.rechunk_for_scaling)")
    lines.append(format_command(["ncatted", "-O", "-h", "-a", "chunk_strategy,global,o,c,scaling", str(out)], args.shell))

    lines.append("\n# 5. verify (_ChunkSizes and time length), then remove intermediates")
    lines.append(format_command(["ncdump", "-hs", str(out)], args.shell))
    lines.append(f"# expect: time = {nt} ; and _ChunkSizes = {', '.join(str(v) for v in final.values())}")
    lines.append(remove_command(intermediates, args.shell))
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Print the NCO/nccopy commands to merge globsim interpolated files from "
                    "disjoint time periods, chunked for single-station reads.")
    parser.add_argument("inputs", nargs="+", type=Path,
                        help="2+ interpolated files, or 2+ interpolation output directories (one per time period)")
    parser.add_argument("-o", "--output-dir", type=Path, required=True, help="directory for the merged files")
    parser.add_argument("--tmp-dir", type=Path, default=None,
                        help="directory for intermediate files (default: OUTPUT_DIR/merge_tmp)")
    parser.add_argument("--kinds", nargs="+", choices=KIND_ORDER, default=["sa", "sf", "pl", "pl_sur"],
                        help="which file kinds to merge when inputs are directories (default: sa sf pl pl_sur)")
    parser.add_argument("--target-mb", type=float, default=SCL_CHK_BYTES / 1024**2,
                        help="target size of one output chunk in MB (default: %(default)s)")
    parser.add_argument("--memory-gb", type=float, default=8,
                        help="memory nccopy may use for its copy buffer + chunk cache (default: %(default)s)")
    parser.add_argument("--deflate", type=int, default=1, choices=range(0, 10),
                        help="deflate level for the final file, 0 = none (default: %(default)s)")
    parser.add_argument("--shell", choices=("bash", "powershell"), default="bash",
                        help="quoting style for the printed commands (default: %(default)s)")
    args = parser.parse_args(argv)

    if len(args.inputs) < 2:
        parser.error("provide at least 2 files or 2 directories")

    tmp_dir = args.tmp_dir or args.output_dir / "merge_tmp"
    groups = find_input_files(args.inputs, args.kinds if all(p.is_dir() for p in args.inputs) else None)
    if not groups:
        parser.error("no globsim interpolated files found")

    mkdir = (["New-Item", "-ItemType", "Directory", "-Force"] if args.shell == "powershell" else ["mkdir", "-p"])
    print(format_command(mkdir + [str(args.output_dir), str(tmp_dir)], args.shell))
    for key, paths in groups.items():
        print()
        print("\n".join(plan_group(key, paths, args.output_dir, tmp_dir, args)))


if __name__ == "__main__":
    main()
