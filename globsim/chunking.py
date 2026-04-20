import netCDF4 as nc
import subprocess
import os

from globsim.memsafe import require_memory_overhead

MIN_CHK_BYTES = 1 * 1024 * 1024  # 1MB 
INTP_CHK_BYTES = 20 * 1024 * 1024  # 20MB 
SCL_CHK_BYTES = 10 * 1024 * 1024  # 10MB 


def _already_chunked_for_scaling(ncf:nc.Dataset | str, bytes_per_element:int=4) -> bool:
    """Check if the file is chunked for scaling (all levels, one station, max time) """
    min_elements = SCL_CHK_BYTES // bytes_per_element
    opened_here = False
    
    if isinstance(ncf, str):
        ncf = nc.Dataset(ncf)
        opened_here = True

    try:
        chunking = None
        chk_dims = None
        
        for _, var in ncf.variables.items():
            chk = var.chunking()
            if chk == 'contiguous':
                continue
            else:
                chk_dims = var.dimensions
                if 'time' not in chk_dims or 'station' not in chk_dims:
                    continue
                chunking = chk
                break

        if chunking is None:
            return False
        
        for dim, cksz in zip(chk_dims, chunking):
            if dim == 'station' and cksz != 1:
                return False
            elif dim == 'level':
                if cksz < ncf.dimensions['level'].size:
                    return False 
                else:
                    min_elements //= cksz
            #elif dim == 'time' and cksz != ncf.dimensions['time'].size and cksz < min_elements:
            #    return False
            
        return True

    finally:
        if opened_here:
            ncf.close()


def find_chunking(ncf:nc.Dataset) -> dict | bool:
    """Find the chunking of the first variable in the file that is not contiguous."""
    opened_here = False
    
    if isinstance(ncf, str):
        ncf = nc.Dataset(ncf)
        opened_here = True
    try:
        chunking = None
        chk_dims = None
        
        for _, var in ncf.variables.items():
            chk = var.chunking()
            if chk == 'contiguous':
                continue
            else:
                chk_dims = var.dimensions
                if 'time' not in chk_dims or 'station' not in chk_dims:
                    continue
                chunking = chk
                break

        if chunking is None:
            return False
        else:
            return dict(zip(chk_dims, chunking))
        
    finally:
        if opened_here:
            ncf.close()


def calculate_chunks_for_interpolation_writing(nt, nl, ns, bytes_per_element=4) -> tuple[int, int, int]:
    """ Calculate chunk sizes for netCDF variables based on the number of time steps (nt), levels (nl), and stations (ns).
    The goal is to optimize read/write performance by keeping chunk sizes reasonable. """
    target_elements = INTP_CHK_BYTES // bytes_per_element  # target number of elements per chunk
    
    c_lev = max(nl, 1)  # if there are no levels, we set c_lev to 1 to avoid division by zero
    c_stat = min(ns, 16)
    
    c_time = int(min(nt, target_elements // (c_lev * c_stat)))
    
    return (c_time, c_lev, c_stat)


def calculate_chunks_for_scaling(nt, nl, ns, bytes_per_element=4) -> tuple[int, int, int]:
    """ Calculate chunk sizes for netCDF variables based on the number of time steps (nt), levels (nl), and stations (ns).
    For scaling, we want all levels in one chunk, one station per chunk, and the rest in time. """
    nl = max(nl, 1)  # if there are no levels, we set nl to 1 to avoid division by zero
    c_stat = 1  # one station per chunk
    c_lev = nl  # all levels in one chunk
    target_elements = SCL_CHK_BYTES // bytes_per_element  # target number of elements per chunk
    
    # Time chunk size is whatever is left after accounting for levels and station
    c_time = int(min(nt, target_elements // (c_lev * c_stat)))
    
    return (c_time, c_lev, c_stat)


def rechunk_for_scaling(f: str):
    """Rechunk a file if it is not already chunked for scaling."""
    if _already_chunked_for_scaling(f):
        print(f"{f} is already chunked for scaling. Skipping rechunking.")
        return
    
    # rename, chunk, verify, delete old file
    # rename input file to a temporary name to avoid overwriting if something goes wrong
    temp_file = f.replace(".nc", "_original.nc")
    output_file = f.replace(".nc", "_chk.nc")
    os.rename(f, temp_file)

    require_memory_overhead(min_gb=2*2)  # Ensure we have enough memory overhead before starting rechunking
    _rechunk_netcdf_for_scaling(temp_file, output_file, memory_buffer_gb=2)
    
    # verify new file reads correctly and has the expected chunking
    try:
        with nc.Dataset(output_file) as ncf:
            if not _already_chunked_for_scaling(ncf):
                print(f"Error: Rechunked file {output_file} is not chunked for scaling as expected.")
                os.rename(temp_file, f)  
            else:
                print(f"Successfully rechunked {f} for scaling.")
                os.remove(temp_file)  
                os.rename(output_file, f) 
    except Exception as e:
        print(f"Error reading rechunked file {output_file}: {e}")
        os.rename(temp_file, f)  
        return


def _rechunk_netcdf_for_scaling(input_file:str, output_file:str, memory_buffer_gb:int=1):
    """
    Rechunks a NetCDF4 file using the nccopy utility.
    
    Args:
        input_file (str): Path to the source file.
        output_file (str): Path where the rechunked file will be saved.
        chunk_spec (str): Chunking definition (e.g., "var1/1,10,10;var2/1,10,10").
        memory_buffer_gb (int): Size of the copy buffer in Gigabytes.
    """
    with nc.Dataset(input_file) as ncf:
        nt = ncf.dimensions['time'].size
        nl = ncf.dimensions['level'].size if 'level' in ncf.dimensions else 0
        ns = ncf.dimensions['station'].size
        print(f"Input file dimensions: time={nt}, level={nl}, station={ns}")
    
    chunks = find_chunking(input_file)
    c_lev = max(1, nl)
    minimum_time = MIN_CHK_BYTES // (4 * c_lev)  
    if chunks is not False: 
        print(f"Existing chunking found: {chunks}")
        existing_time = chunks.get('time', minimum_time)  # match time slice so nccopy doesn't run out of mem
        c_time = max(existing_time, minimum_time)  # Avoid tiny chunks if file is poorly chunked
        c_time = min(nt, c_time)  # Don't create chunks larger than the dimension size
        c_stat = 1
    else:
        c_time, c_lev, c_stat = calculate_chunks_for_scaling(nt, nl, ns)

    if nl == 0:
        chunk_spec = f"time/{c_time},station/{c_stat}"
    else:
        chunk_spec = f"time/{c_time},level/{c_lev},station/{c_stat}"
    # Command: nccopy -c [chunks] -m [buffer] input output
    # -c defines the new chunk shapes
    # -m defines the size of the copy buffer (helps with large files)
    cmd = [
        "nccopy", 
        "-k", "netCDF-4", # Ensure output is NetCDF4
        "-c", chunk_spec, 
        "-d", "1",
        "-m", f"{memory_buffer_gb}G", 
        input_file, 
        output_file
    ]
    print(" ".join(cmd))
    print(f"Starting rechunking: {input_file} -> {output_file}")
    success = None
    try:
        subprocess.run(cmd, check=True)
        print("Rechunking complete.")
        success = True
    except subprocess.CalledProcessError as e:
        print(f"Error during rechunking: {e}")
        success = False

    if success:
        with nc.Dataset(output_file, 'a') as ncf:
            ncf.setncattr('chunk_strategy', 'scaling')
        
    return success


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python chunking.py <file.nc>")
    else:
        rechunk_for_scaling(sys.argv[1])