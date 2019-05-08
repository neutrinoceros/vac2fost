"""Test input and output facilities when current dir is either the inputidir (origin) or elsewhere..."""

import os
import pytest
from pathlib import Path
from amrvac_pywrap import secure_chdir
from vac2fost import Interface

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
        Interface(config_file=test_dir/"sample/vac2fost_conf_nonaxisym.nml",
                  mcfost_verbose=True)
