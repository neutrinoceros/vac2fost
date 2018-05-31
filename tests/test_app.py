'''Test basic consistency of main() called through python.'''

import os
import pathlib
import subprocess
import pytest

import f90nml
import astropy.io.fits as pyfits

from vac2fost import main as app

here = pathlib.Path(__file__).absolute().parent


@pytest.mark.incremental #each test is run only if the previous one passed
class TestPyScripting():
    def test_python_call(self):
        app(str(here/'sample/vac2fost_conf.nml'))

    def test_format(self):
        f = pyfits.open(here.parent / 'hd142527_dusty0000.fits')[0]
        opt = f90nml.read(here / 'sample/vac2fost_conf.nml')['mcfost_list']
        assert f.data.shape[1:] == (opt['nphi'], 2*opt['nz']+1, opt['nr'])
