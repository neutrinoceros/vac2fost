'''Test basic consistency of main() called through python.'''
import pathlib
import shutil
import pytest

import f90nml
import astropy.io.fits as pyfits

from vac2fost import main as app

testdir = pathlib.Path(__file__).absolute().parent.parent

@pytest.mark.incremental #each test is run only if the previous one passed
class TestPyScripting():
    output_dir = testdir / "output/TestPyScripting"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    def test_python_call(self):
        app(
            str(testdir/"sample/vac2fost_conf.nml"),
            output_dir=__class__.output_dir,
            mcfost_verbose=True,
        )

    def test_format(self):
        f = pyfits.open(__class__.output_dir / 'hd142527_dusty0000.fits')[0]
        opt = f90nml.read(testdir / 'sample/vac2fost_conf.nml')['mcfost_output']
        assert f.data.shape[1:] == (opt['n_az'], opt['nz'], opt['n_rad'])

    def test_python_call_multiple(self):
        app(
            str(testdir / "sample/vac2fost_conf.nml"),
            output_dir=__class__.output_dir,
            mcfost_verbose=True,
            nums=[0, 1, 2]
        )
        for n in (0, 1, 2):
            assert (__class__.output_dir / f"hd142527_dusty{str(n).zfill(4)}.fits").exists()
