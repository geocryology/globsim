"""
Central registry of reanalysis backends.
Adding a new reanalysis = adding one entry here + writing the three classes.
"""
from dataclasses import dataclass, field
from importlib import import_module
from typing import Type, Optional


def _resolve_class(dotted_path: str):
    """Import 'globsim.scale.ERA5scale.ERA5scale' → <class ERA5scale>"""
    module_path, class_name = dotted_path.rsplit(".", 1)
    module = import_module(module_path)
    return getattr(module, class_name)


def parse_enabled(args_d: "list[str] | None") -> dict[str, bool]:
    """Convert CLI's -d list into {backend_key: bool} dict.
    
    If args_d is None, nothing is enabled (enabled_by_default=False).
    If args_d contains backend keys, only those are enabled.
    """
    if args_d is None:
        return {}
    return {key: (key in args_d) for key in BACKENDS}


@dataclass
class ReanalysisBackend:
    """One reanalysis product (e.g. ERA5, MERRA-2)."""
    name: str                          # Human-readable, e.g. "ERA-5"
    key: str                           # Short key used in CLI / config, e.g. "ERA5"
    download_cls: str                  # Dotted import path, e.g. "globsim.download.ERA5MonthlyDownload"
    interpolate_cls: str               # e.g. "globsim.interpolate.ERA5interpolate"
    scale_cls: str                     # e.g. "globsim.scale.ERA5scale"
    download_args: dict = field(default_factory=dict)  # extra kwargs for download 
    enabled_by_default: bool = False # whether to run this backend if no -d option is given


BACKENDS: dict[str, ReanalysisBackend] = {}

def register(backend: ReanalysisBackend):
    BACKENDS[backend.key] = backend

def get_backend(key: str) -> ReanalysisBackend:
    return BACKENDS[key]

def all_keys() -> list[str]:
    return list(BACKENDS.keys())

# ── Register existing backends ────────────────────────────────────
register(ReanalysisBackend(
    name="ERA-5",
    key="ERA5",
    download_cls="globsim.download.era5_monthly.ERA5MonthlyDownload",
    interpolate_cls="globsim.interpolate.ERA5interpolate.ERA5interpolate",
    scale_cls="globsim.scale.ERA5scale.ERA5scale",
    download_args={"ens": False},
))

register(ReanalysisBackend(
    name="ERA-5 Ensemble",
    key="ERA5ENS",
    download_cls="globsim.download.ERA5Edownload.ERA5Edownload",
    interpolate_cls="globsim.interpolate.ERA5interpolate.ERA5EnsembleInterpolate",
    scale_cls="globsim.scale.ERA5Escale.ERA5Escale",
))

register(ReanalysisBackend(
    name="MERRA-2",
    key="MERRA",
    download_cls="globsim.download.MERRAdownload.MERRAdownload",
    interpolate_cls="globsim.interpolate.MERRAinterpolate.MERRAinterpolate",
    scale_cls="globsim.scale.MERRAscale.MERRAscale",
))

register(ReanalysisBackend(
    name="JRA-55",
    key="JRA",
    download_cls="globsim.download.JRAdownload.JRAdownload",
    interpolate_cls="globsim.interpolate.JRAinterpolate.JRAinterpolate",
    scale_cls="globsim.scale.JRA3Qscale.JRA55",
))

register(ReanalysisBackend(
    name="JRA-3Q",
    key="JRA3Q",
    download_cls="globsim.download.JRA3Qdownload.J3QD",
    interpolate_cls="globsim.interpolate.JRA3Qinterpolate.J3QI",
    scale_cls="globsim.scale.JRA3Qscale.J3QS",
))

register(ReanalysisBackend(
    name="JRA-3Q Gaussian",
    key="JRA3QG",
    download_cls="globsim.download.JRA3Qdownload.J3QgD",
    interpolate_cls="globsim.interpolate.JRA3Qinterpolate.J3QgI",
    scale_cls="globsim.scale.JRA3Qscale.J3QgS",
))