import pickle
import pathlib
import os
import numpy as np
from astropy.io import fits as pyfits

from vac2fost import main as app

here = pathlib.Path(__file__).absolute().parent
outdir = here / 'output/test_regression'
gridfile = outdir / 'mcfost_grid.fits.gz'
if gridfile.exists():
    os.remove(gridfile)

out = app(
    str(here/'sample/vac2fost_conf.nml'),
    output_dir=outdir
)

class TestRegression:
    def test_mcfost_conf(self):
        with open(here/'ref/mcfost_conf.para') as fi:
            ref_lines = fi.readlines()
        with open(outdir/'mcfost_conf.para') as fi:
            new_lines = fi.readlines()
        for n, r in zip(new_lines[:-2], ref_lines):
            assert n==r

    def test_target_grid(self):
        ref = pyfits.open(here/'ref/mcfost_grid.fits.gz')[0].data
        new = pyfits.open(outdir/'mcfost_grid.fits.gz')[0].data
        assert (ref == new).all()

    def test_image(self):
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        data = pyfits.open(out['finame'])[0].data[0]

        ref = pyfits.open(here/'ref/hd142527_dusty0000.fits')[0].data[0]
        diff = data - ref
        assert np.abs(diff).max() < 1e-15

    def test_out(self):
        out_ref = pickle.load(open(here/'ref/main_out.p', 'rb'))
        for k,v in out.items():
            if isinstance(v, np.ndarray):
                assert (v == out_ref[k]).all()
