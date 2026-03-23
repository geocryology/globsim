from multiprocessing.dummy import Pool as ThreadPool
from globsim.registry import BACKENDS, _resolve_class


def _selected_backends(enabled: dict[str, bool]) -> list:
    """Return backend entries where enabled[key] is True."""
    return [b for b in BACKENDS.values()
            if enabled.get(b.key, b.enabled_by_default)]


def GlobsimDownload(pfile, multithread=True, **enabled):
    """Download data from selected reanalyses."""
    objects = []
    for b in _selected_backends(enabled):
        cls = _resolve_class(b.download_cls)
        obj = cls(pfile, **b.download_args) if b.download_args else cls(pfile)
        objects.append(obj)

    if multithread and len(objects) > 1:
        pool = ThreadPool(len(objects))
        pool.map(lambda ob: ob.retrieve(), objects)
        pool.close()
        pool.join()
    else:
        for ob in objects:
            ob.retrieve()


def GlobsimInterpolateStation(ifile, **enabled_and_kwargs):
    """Interpolate re-analysis grids to station points."""
    # Separate backend toggles from passthrough kwargs
    enabled = {k: v for k, v in enabled_and_kwargs.items() if k in BACKENDS}
    kwargs  = {k: v for k, v in enabled_and_kwargs.items() if k not in BACKENDS}

    for b in _selected_backends(enabled):
        cls = _resolve_class(b.interpolate_cls)
        cls(ifile, **kwargs).process()


def GlobsimScale(sfile, **enabled):
    """Scale interpolated data to near-surface fluxes."""
    for b in _selected_backends(enabled):
        cls = _resolve_class(b.scale_cls)
        cls(sfile).process()