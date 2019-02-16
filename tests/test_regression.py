import pickle
from pathlib import Path
import os
import shutil
import numpy as np
from astropy.io import fits as pyfits

from vac2fost import main as app
from vac2fost.vac2fost import DETECTED_MCFOST_VERSION

x,y,z = DETECTED_MCFOST_VERSION
assert x == 3 and y == 0
if z < 35:
    REFVER = "3.0.34"
else:
    REFVER = "3.0.35"

test_dir = Path(__file__).parent.resolve()
REFOUT_DIR = test_dir / f"ref/{REFVER}"
OUT = test_dir/"output"

def instanciate_interface(conffile):
    outdir = OUT / f"test_reg_{Path(conffile).stem}"
    if outdir.is_dir():
        shutil.rmtree(outdir)
    itf = app(test_dir/"sample"/conffile, output_dir=outdir)
    return itf

# to regold tests
save_keys = ['sim_conf',
             'input_grid', 'output_grid',
             'new_2D_arrays', 'new_3D_arrays',
             'dust_binning_mode'
]

def regold(itf, reffile):
    with open(reffile, mode="wb") as file:
        out = {k: itf.__getattribute__(k) for k in save_keys}
        pickle.dump(out, file)

class TestRegressionMain:
    subrefdir = REFOUT_DIR / "default"
    itf = instanciate_interface(conffile="vac2fost_conf.nml")
    itf.tag = itf._base_args['config_file'].stem

    def test_mcfost_conf(self):
        itf = __class__.itf
        with open(__class__.subrefdir / "mcfost_conf.para") as fi:
            ref_lines = fi.readlines()
        with open(itf.io["out"].directory / "mcfost_conf.para") as fi:
            new_lines = fi.readlines()
        for n, r in zip(new_lines[:-2], ref_lines):
            assert n == r

    def test_target_grid(self):
        itf = __class__.itf
        ref = pyfits.open(__class__.subrefdir / "mcfost_grid.fits.gz")[0].data
        new = pyfits.open(itf.io["out"].directory/'mcfost_grid.fits.gz')[0].data
        np.testing.assert_array_equal(ref, new)

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        #regold(itf, reffile)

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf.dust_binning_mode == out_ref['dust_binning_mode']
        assert itf.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_array_equal(itf.output_grid['zg'], out_ref['output_grid']['zg'])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)

    def test_image(self):
        itf = __class__.itf
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        fipath = itf.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]
        ref = pyfits.open(__class__.subrefdir / "hd142527_dusty0000.fits")[0].data[0]
        np.testing.assert_array_equal(data, ref)

class TestRegressionNonAxisym:
    subrefdir = REFOUT_DIR / "nonaxisym"
    itf = instanciate_interface(conffile="vac2fost_conf_nonaxisym.nml")
    itf.tag = itf._base_args['config_file'].stem

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        #regold(itf, reffile)

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf.dust_binning_mode == out_ref['dust_binning_mode']
        assert itf.sim_conf == out_ref['sim_conf']
        np.testing.assert_array_equal(itf.input_grid['rv'], out_ref['input_grid']['rv'])
        np.testing.assert_array_equal(itf.input_grid['phiv'], out_ref['input_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rv'], out_ref['output_grid']['rv'])
        np.testing.assert_array_equal(itf.output_grid['phiv'], out_ref['output_grid']['phiv'])
        np.testing.assert_array_equal(itf.output_grid['rg'], out_ref['output_grid']['rg'])
        np.testing.assert_array_equal(itf.output_grid['phig'], out_ref['output_grid']['phig'])
        np.testing.assert_array_equal(itf.output_grid['zg'], out_ref['output_grid']['zg'])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref['new_2D_arrays'], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref['new_3D_arrays'], rtol=1e-15)

    def test_image(self):
        itf = __class__.itf
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        fipath = itf.io['out'].filepath
        data = pyfits.open(fipath)[0].data[0]
        ref = pyfits.open(__class__.subrefdir / "hd142527_dusty0000_nonaxisym.fits")[0].data[0]
        np.testing.assert_array_equal(data, ref)

class TestRegressionAutoGasOnly:
    subrefdir = REFOUT_DIR / "autogasonly"
    itf = instanciate_interface(conffile="autogasonly/rwi.nml")
    itf.tag = itf._base_args['config_file'].stem

    def test_mcfost_conf(self):
        itf = __class__.itf
        with open(__class__.subrefdir / "mcfost_conf.para") as fi:
            ref_lines = fi.readlines()
        with open(itf.io["out"].directory/"mcfost_conf.para") as fi:
            new_lines = fi.readlines()
        for n, r in zip(new_lines[:-2], ref_lines):
            assert n == r

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        #regold(itf, reffile)

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf.dust_binning_mode == out_ref["dust_binning_mode"]
        assert itf.sim_conf == out_ref["sim_conf"]
        np.testing.assert_array_equal(itf.input_grid["rv"], out_ref["input_grid"]["rv"])
        np.testing.assert_array_equal(itf.input_grid["phiv"], out_ref["input_grid"]["phiv"])
        np.testing.assert_array_equal(itf.output_grid["rv"], out_ref["output_grid"]["rv"])
        np.testing.assert_array_equal(itf.output_grid["phiv"], out_ref["output_grid"]["phiv"])
        np.testing.assert_array_equal(itf.output_grid["rg"], out_ref["output_grid"]["rg"])
        np.testing.assert_array_equal(itf.output_grid["phig"], out_ref["output_grid"]["phig"])
        np.testing.assert_array_equal(itf.output_grid["zg"], out_ref["output_grid"]["zg"])
        np.testing.assert_allclose(itf.new_2D_arrays, out_ref["new_2D_arrays"], rtol=1e-25)
        np.testing.assert_allclose(itf.new_3D_arrays, out_ref["new_3D_arrays"], rtol=1e-15)
