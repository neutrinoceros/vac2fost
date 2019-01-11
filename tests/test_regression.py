import pickle
import pathlib
import os
import shutil
import numpy as np
from astropy.io import fits as pyfits

from vac2fost import main as app

test_dir = pathlib.Path(__file__).absolute().parent
outdir = test_dir / 'output/test_regression'
gridfile = outdir / 'mcfost_grid.fits.gz'
if outdir.is_dir():
    shutil.rmtree(outdir)

itf = app(test_dir/'sample/vac2fost_conf.nml', output_dir=outdir)

class TestRegression:
    def test_mcfost_conf(self):
        with open(test_dir/'ref/mcfost_conf.para') as fi:
            ref_lines = fi.readlines()
        with open(outdir/'mcfost_conf.para') as fi:
            new_lines = fi.readlines()
        for n, r in zip(new_lines[:-2], ref_lines):
            assert n == r

    def test_target_grid(self):
        ref = pyfits.open(test_dir/'ref/mcfost_grid.fits.gz')[0].data
        new = pyfits.open(outdir/'mcfost_grid.fits.gz')[0].data
        np.testing.assert_array_equal(ref, new)

    def test_image(self):
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        fipath = itf.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]

        ref = pyfits.open(test_dir/'ref/hd142527_dusty0000.fits')[0].data[0]
        np.testing.assert_allclose(data, ref, rtol=1e-15)

    def test_out(self):
        out_ref = pickle.load(open(test_dir/'ref/main_out.p', 'rb'))

        #use this to regold the reference file
        # with open(test_dir/'ref/main_out.p', 'wb') as file:
        #     save_keys = ['sim_conf',
        #                 'input_grid', 'output_grid',
        #                  'new_2D_arrays', 'new_3D_arrays',
        #                  '_dbm'
        # ]
        #     out = {k: itf.__getattribute__(k) for k in save_keys}
        #     pickle.dump(out, file)

        assert itf.dust_binning_mode == out_ref['_dbm']
        assert itf.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)

    def test_out_non_axisym(self):
        outdir = test_dir / 'output/test_regression_non_axisym'
        if outdir.is_dir():
            shutil.rmtree(outdir)
        itf = app(test_dir/'sample/vac2fost_conf_nonaxisym.nml', output_dir=outdir)
        out_ref = pickle.load(open(test_dir/'ref/main_out_non_axisym.p', 'rb'))

        #use this to regold the reference file
        # with open(test_dir/'ref/main_out_non_axisym.p', 'wb') as file:
        #     save_keys = ['sim_conf',
        #                  'input_grid', 'output_grid',
        #                  'new_2D_arrays', 'new_3D_arrays',
        #                  '_dbm'
        #     ]
        #     out = {k: itf.__getattribute__(k) for k in save_keys}
        #     pickle.dump(out, file)

        assert itf.dust_binning_mode == out_ref['_dbm']
        assert itf.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)
