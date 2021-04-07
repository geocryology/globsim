""" Save us the trouble of importing ALL the ESMF machinery on scripts that load interpolating functions """
import importlib


class LazyLoader:
    def __init__(self, lib_name):
        self.lib_name = lib_name
        self._mod = None

    def __getattr__(self, name):
        if self._mod is None:
            self._mod = importlib.import_module(self.lib_name)
            
        return getattr(self._mod, name)