"""Instantiate a module-wise logger"""
import logging

v2flogger = logging.getLogger("vac2fost")
handle = logging.StreamHandler()
handle.setFormatter(logging.Formatter("%(name)s | %(levelname)s - %(message)s"))
v2flogger.addHandler(handle)
v2flogger.propagate = False
