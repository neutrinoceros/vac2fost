"""Test basic consistency of main() called through python."""
import pathlib
import shutil

import f90nml
import astropy.io.fits as pyfits

from vac2fost import main as app
from vac2fost.logger import v2flogger as log

testdir = pathlib.Path(__file__).absolute().parent.parent


class TestPyScripting:
    output_dir = testdir / "output/TestPyScripting"
    if output_dir.exists():
        shutil.rmtree(output_dir)

    def test_python_call(self):
        log.setLevel(10)
        app(str(testdir / "sample/vac2fost_conf.nml"), output_dir=__class__.output_dir)

    def test_format(self):
        f = pyfits.open(__class__.output_dir / "hd142527_dusty0000.fits")[0]
        opt = f90nml.read(testdir / "sample/vac2fost_conf.nml")["mcfost_output"]
        assert f.data.shape[1:] == (opt["n_az"], opt["nz"], opt["n_rad"])

    def test_python_call_multiple(self):
        log.setLevel(10)
        app(
            str(testdir / "sample/vac2fost_conf.nml"),
            output_dir=__class__.output_dir,
            override={"amrvac_input": {"nums": [0, 1, 2]}},
        )
        for n in (0, 1, 2):
            assert (__class__.output_dir / f"hd142527_dusty{str(n).zfill(4)}.fits").exists()
