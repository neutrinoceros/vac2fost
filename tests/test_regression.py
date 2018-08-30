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

itf = app(
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
        itf.write_output()
        fipath = itf.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]

        ref = pyfits.open(here/'ref/hd142527_dusty0000.fits')[0].data[0]
        diff = data - ref
        assert np.abs(diff).max() < 1e-15

    def test_out(self):
        out_ref = pickle.load(open(here/'ref/main_out.p', 'rb'))
        save_keys = ['sim_conf',
                     'input_grid', 'output_grid',
                     'new_2D_arrays', 'new_3D_arrays',
                     '_dbm'
        ]
        out = {k: itf.__getattribute__(k) for k in save_keys}
        # use this to regold the reference file
        #with open(here/'ref/main_out2.p', 'wb') as fo:
        #    pickle.dump(out, fo)

        print(f"\n\nDBM, new vs ref: {out['_dbm']}, {out_ref['_dbm']}")
        assert deep_equality(out, out_ref)

def deep_equality(a, b) -> bool:
    res = False
    if isinstance(a, dict):
        for k in a.keys():
            if not deep_equality(a[k], b[k]): break
        else: res = True
    elif isinstance(a, list):
        for a1,b1 in zip(a,b):
            if not deep_equality(a1, b1): break
        else: res = True
    elif isinstance(a, np.ndarray):
        res = (a == b).all()
    else:
        try:
            res = (a == b)
        except ValueError:
            raise ValueError(f'raised for object type {type(a)}')

    return res

