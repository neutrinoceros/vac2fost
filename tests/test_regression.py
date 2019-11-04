import pickle
from pathlib import Path
import shutil
import numpy as np
from astropy.io import fits

from vac2fost import main as app
from vac2fost.logger import v2flogger
import pytest

test_dir = Path(__file__).parent.resolve()
REFOUT_DIR = test_dir / "ref"
OUT = test_dir/"output"

def instanciate_interface(conffile, **kwargs):
    outdir = OUT / f"test_reg_{Path(conffile).stem}"
    if outdir.is_dir():
        shutil.rmtree(outdir)
    override = {"flags": {k: v for k,v in kwargs.items()}}
    itf = app(test_dir/"sample"/conffile, override, output_dir=outdir)
    return itf

# to regold tests
def regold(itf, reffile, morekeys:list = None):
    save_keys = [
        "amrvac_conf",
        "input_grid", "output_grid",
        "_dust_binning_mode"
    ]
    if morekeys is not None:
        save_keys += morekeys
    with open(reffile, mode="wb") as file:
        out = {k: itf.__getattribute__(k) for k in save_keys}
        out.update({"output_ndarray": itf.get_output_ndarray()})
        pickle.dump(out, file)

class TestRegressionMain:
    subrefdir = REFOUT_DIR / "default"
    v2flogger.setLevel(10) # debug
    itf = instanciate_interface(conffile="vac2fost_conf.nml", read_gas_velocity=True)
    itf.preroll_mcfost()
    itf.tag = itf.conf_file.stem

    def test_mcfost_conf(self):
        itf = __class__.itf
        with open(itf.io.OUT.directory / "mcfost_conf.para") as fi:
            new_lines = fi.readlines()
        #with open(__class__.subrefdir / "mcfost_conf.para", mode="wt") as fi: #regold...
        #    fi.write("".join(new_lines))
        with open(__class__.subrefdir / "mcfost_conf.para") as fi:
            ref_lines = fi.readlines()
        for n, r in zip(new_lines[2:-3], ref_lines[2:]):
            assert n == r

    def test_target_grid(self):
        itf = __class__.itf
        ref_file = __class__.subrefdir / "mcfost_grid.fits.gz"
        new_file = itf.io.OUT.directory/"mcfost_grid.fits.gz"
        ref = fits.open(ref_file)[0].data
        new = fits.open(new_file)[0].data
        #shutil.copyfile(new_file, ref_file) # to regold ...
        np.testing.assert_allclose(ref, new, rtol=1e-15)

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        #regold(itf, reffile, morekeys=["new_3D_gas_velocity"])

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf._dust_binning_mode == out_ref["_dust_binning_mode"]
        assert itf.amrvac_conf == out_ref["amrvac_conf"]
        np.testing.assert_array_equal(itf.input_grid["ticks_r"], out_ref["input_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.input_grid["ticks_phi"], out_ref["input_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["ticks_r"], out_ref["output_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.output_grid["ticks_phi"], out_ref["output_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_r"], out_ref["output_grid"]["z-slice_r"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_phi"], out_ref["output_grid"]["z-slice_phi"])
        np.testing.assert_allclose(itf.output_grid["phi-slice_z"], out_ref["output_grid"]["phi-slice_z"], rtol=1e-15)
        np.testing.assert_allclose(itf.get_output_ndarray(), out_ref["output_ndarray"], rtol=5e-14)
        np.testing.assert_allclose(itf.new_3D_gas_velocity, out_ref["new_3D_gas_velocity"], rtol=5e-14)

    def test_image(self):
        itf = __class__.itf
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        ref_file = __class__.subrefdir / "hd142527_dusty0000.fits"
        new_file = itf.io.OUT.filepath
        #shutil.copyfile(new_file, ref_file) # to regold ...
        ref = fits.open(ref_file)[0].data[0]
        new = fits.open(new_file)[0].data[0]
        np.testing.assert_allclose(new, ref, rtol=5e-14)


@pytest.mark.skip(reason="wip")
class TestRegressionMutliNums:
    subrefdir = REFOUT_DIR / "multinums"
    itf = instanciate_interface(conffile="vac2fost_conf_quick.nml", nums=[0, 1, 2])
    itf.preroll_mcfost()
    itf.tag = itf.conf_file.stem + "multinums"

    def test_multinums_output(self):
        filename = __class__.itf.io.OUT.filepath.stem[:-4]
        for n in (0, 1, 2):
            new_file = __class__.itf.io.OUT.directory / f"{filename}{str(n).zfill(4)}.fits"
            ref_file = __class__.subrefdir / f"hd142527_dusty{str(n).zfill(4)}.fits"
            assert new_file.exists()
            #shutil.copyfile(new_file, ref_file) # to regold ...
            np.testing.assert_allclose(fits.open(new_file)[0].data, fits.open(ref_file)[0].data, rtol=5e-14)


class TestRegressionNonAxisym:
    subrefdir = REFOUT_DIR / "nonaxisym"
    itf = instanciate_interface(conffile="vac2fost_conf_nonaxisym.nml")
    itf.preroll_mcfost()
    itf.tag = itf.conf_file.stem

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        # nb: the present test forces the computatio of the 3D gas velocity field
        # BY testing or regolding it, despite the fact it's not asked at instanciation
        #regold(itf, reffile, morekeys=["new_3D_gas_velocity"])

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf._dust_binning_mode == out_ref["_dust_binning_mode"]
        assert itf.amrvac_conf == out_ref["amrvac_conf"]
        np.testing.assert_array_equal(itf.input_grid["ticks_r"], out_ref["input_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.input_grid["ticks_phi"], out_ref["input_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["ticks_r"], out_ref["output_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.output_grid["ticks_phi"], out_ref["output_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_r"], out_ref["output_grid"]["z-slice_r"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_phi"], out_ref["output_grid"]["z-slice_phi"])
        np.testing.assert_allclose(itf.output_grid["phi-slice_z"], out_ref["output_grid"]["phi-slice_z"], rtol=1e-15)
        np.testing.assert_allclose(itf.get_output_ndarray(), out_ref["output_ndarray"], rtol=5e-14)
        np.testing.assert_allclose(itf.new_3D_gas_velocity, out_ref["new_3D_gas_velocity"], rtol=5e-14)

    def test_image(self):
        itf = __class__.itf
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        ref_file = __class__.subrefdir / "hd142527_rphi0020.fits"
        new_file = itf.io.OUT.filepath
        ref = fits.open(ref_file)[0].data[0]
        new = fits.open(new_file)[0].data[0]
        #shutil.copyfile(new_file, ref_file) # to regold ...
        np.testing.assert_allclose(new, ref, rtol=5e-14)


class TestRegressionAutoGasOnly:
    subrefdir = REFOUT_DIR / "autogasonly"
    itf = instanciate_interface(conffile="autogasonly/rwi.nml")
    itf.preroll_mcfost()
    itf.tag = itf.conf_file.stem

    def test_mcfost_conf(self):
        itf = __class__.itf
        with open(itf.io.OUT.directory / "mcfost_conf.para") as fi:
            new_lines = fi.readlines()
        #with open(__class__.subrefdir / "mcfost_conf.para", mode="wt") as fi: #regold...
        #    fi.write("".join(new_lines))
        with open(__class__.subrefdir / "mcfost_conf.para") as fi:
            ref_lines = fi.readlines()
        for n, r in zip(new_lines[2:-3], ref_lines[2:]):
            assert n == r

    def test_out(self):
        itf = __class__.itf
        reffile = __class__.subrefdir / f"{itf.tag}.p"
        #regold(itf, reffile)

        out_ref = pickle.load(open(reffile, mode="rb"))
        assert itf._dust_binning_mode == out_ref["_dust_binning_mode"]
        assert itf.amrvac_conf == out_ref["amrvac_conf"]
        np.testing.assert_array_equal(itf.input_grid["ticks_r"], out_ref["input_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.input_grid["ticks_phi"], out_ref["input_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["ticks_r"], out_ref["output_grid"]["ticks_r"])
        np.testing.assert_array_equal(itf.output_grid["ticks_phi"], out_ref["output_grid"]["ticks_phi"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_r"], out_ref["output_grid"]["z-slice_r"])
        np.testing.assert_array_equal(itf.output_grid["z-slice_phi"], out_ref["output_grid"]["z-slice_phi"])
        np.testing.assert_array_equal(itf.output_grid["phi-slice_z"], out_ref["output_grid"]["phi-slice_z"])
        np.testing.assert_allclose(itf.get_output_ndarray(), out_ref["output_ndarray"], rtol=1e-15)
