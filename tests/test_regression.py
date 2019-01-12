import pickle
import pathlib
import os
import shutil
import numpy as np
from astropy.io import fits as pyfits

from vac2fost import main as app

test_dir = pathlib.Path(__file__).absolute().parent
outdir = test_dir / 'output/test_regression'
if outdir.is_dir():
    shutil.rmtree(outdir)
itf = app(test_dir/'sample/vac2fost_conf.nml', output_dir=outdir)

outdir2 = test_dir / 'output/test_regression_non_axisym'
if outdir2.is_dir():
    shutil.rmtree(outdir2)
itf2 = app(test_dir/'sample/vac2fost_conf_nonaxisym.nml', output_dir=outdir2)

# to regold tests
save_keys = ['sim_conf',
             'input_grid', 'output_grid',
             'new_2D_arrays', 'new_3D_arrays',
             'dust_binning_mode'
]

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

    def test_out(self):
        #use this to regold the reference file
        # with open(test_dir/'ref/main_out.p', 'wb') as file:
        #     out = {k: itf.__getattribute__(k) for k in save_keys}
        #     pickle.dump(out, file)

        out_ref = pickle.load(open(test_dir/'ref/main_out.p', 'rb'))
        assert itf.dust_binning_mode == out_ref['dust_binning_mode']
        assert itf.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)

    def test_image(self):
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        fipath = itf.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]
        ref = pyfits.open(test_dir/'ref/hd142527_dusty0000.fits')[0].data[0]
        np.testing.assert_array_equal(data, ref)

    def test_out_non_axisym(self):
        #use this to regold the reference file
        # with open(test_dir/'ref/main_out_non_axisym.p', 'wb') as file:
        #     out = {k: itf.__getattribute__(k) for k in save_keys}
        #     pickle.dump(out, file)

        out_ref = pickle.load(open(test_dir/'ref/main_out_non_axisym.p', 'rb'))
        assert itf2.dust_binning_mode == out_ref['dust_binning_mode']
        assert itf2.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf2.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf2.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf2.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf2.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf2.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf2.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_allclose(itf2.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf2.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)

    def test_image_non_axisym(self):
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf2.write_output()
        fipath = itf2.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]
        ref = pyfits.open(test_dir/'ref/hd142527_dusty0000_nonaxisym.fits')[0].data[0]
        np.testing.assert_array_equal(data, ref)
