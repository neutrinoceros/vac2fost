"""Test input and output facilities when current dir is either the inputidir (origin) or elsewhere..."""

import os
import pytest
from pathlib import Path
from contextlib import contextmanager

from vac2fost.logger import v2flogger as log
from vac2fost import VtuFileInterface as Interface

@contextmanager
def secure_chdir(path):
    '''A context manager that cd back to original location at closing.'''
    origin = os.getcwd()
    os.chdir(path)
    try: yield
    finally: os.chdir(origin)

test_dir = Path(__file__).absolute().parent.parent

sampledir = test_dir / "sample"

class Test_IO:
    def test_input_fromloc(self):
        with secure_chdir(sampledir):
            itf = Interface("vac2fost_conf.nml")
        assert itf.io.IN.directory == Path(sampledir).resolve()

    def test_input_fromhome(self):
        with secure_chdir(os.environ["HOME"]):
            itf = Interface(sampledir/"vac2fost_conf.nml")
        assert itf.io.IN.directory == Path(sampledir).resolve()

    def test_input_fromhome_fakedata(self):
        with pytest.raises(FileNotFoundError):
            with secure_chdir(os.environ["HOME"]):
                Interface("vac2fost_conf.nml")

    def test_output_fromloc(self):
        with secure_chdir(sampledir):
            itf = Interface("vac2fost_conf.nml")
        print(itf.io.OUT.directory)

    def test_path_reading(self):
        """Check that AMRVAC config file can be correctly assessed with a
        relative "origin" argument"""
        log.setLevel(10)
        Interface(test_dir/"sample/vac2fost_conf_nonaxisym.nml")
