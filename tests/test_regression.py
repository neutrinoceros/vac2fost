import pickle
from pathlib import Path
import shutil
import numpy as np
from astropy.io import fits

from vac2fost import main
from vac2fost import VtuFileInterface as Interface
from vac2fost.logger import v2flogger
from conftest import TEST_DATA_DIR, TEST_VTU_DATA_DIR, TEST_DAT_DATA_DIR, TEST_ARTIFACTS_DIR


v2flogger.setLevel(10)




# to regold tests
def regold(itf, reffile, morekeys: list = None):
    save_keys = ["amrvac_conf", "input_grid", "output_grid", "_dust_bin_mode"]
    if morekeys is not None:
        save_keys += morekeys
    out = {k: itf.__getattribute__(k) for k in save_keys}
    out.update({"output_ndarray": itf.get_output_ndarray()})
    if "read_gas_velocity" in itf.conf["flags"]:
        out.update({"output_gas_velocity": itf.get_gas_velocity_ndarray()})
    with open(reffile, mode="wb") as file:
        pickle.dump(out, file)


class AbstractTestRegression:
    def test_out(self):
        itf = self.__class__.itf

        out_ref = pickle.load(open(self.__class__.answerfile, mode="rb"))
        assert itf._dust_bin_mode == out_ref["_dust_bin_mode"]
        assert itf.amrvac_conf == out_ref["amrvac_conf"]
        np.testing.assert_array_equal(itf.input_grid["ticks_r"], out_ref["input_grid"]["ticks_r"])
        np.testing.assert_array_equal(
            itf.input_grid["ticks_phi"], out_ref["input_grid"]["ticks_phi"]
        )
        np.testing.assert_array_equal(itf.output_grid["ticks_r"], out_ref["output_grid"]["ticks_r"])
        np.testing.assert_array_equal(
            itf.output_grid["ticks_phi"], out_ref["output_grid"]["ticks_phi"]
        )
        np.testing.assert_array_equal(
            itf.output_grid["z-slice_r"], out_ref["output_grid"]["z-slice_r"]
        )
        np.testing.assert_array_equal(
            itf.output_grid["z-slice_phi"], out_ref["output_grid"]["z-slice_phi"]
        )
        np.testing.assert_allclose(
            itf.output_grid["phi-slice_z"], out_ref["output_grid"]["phi-slice_z"], rtol=1e-15
        )
        np.testing.assert_allclose(itf.get_output_ndarray(), out_ref["output_ndarray"], rtol=5e-14)
        if "output_gas_velocity" in out_ref.keys():
            np.testing.assert_allclose(
                itf.get_gas_velocity_ndarray(), out_ref["output_gas_velocity"], rtol=5e-14
            )


class AbstractTestImage(AbstractTestRegression):
    def test_image(self):
        itf = self.__class__.itf
        # get the Primary (only image available),
        # and exctract its first 3d array (density field)
        itf.write_output()
        new_file = itf.io.OUT.filepath
        ref_file = self.__class__.answerdir / new_file.name
        # shutil.copyfile(new_file, ref_file) # to regold ...
        ref = fits.open(ref_file)[0].data[0]
        new = fits.open(new_file)[0].data[0]
        np.testing.assert_allclose(new, ref, rtol=5e-14)


# Actual test classes:
class TestRegressionAutoGasOnly(AbstractTestRegression):
    stem = "autogasonly"
    answerdir = TEST_VTU_DATA_DIR / stem
    conffile = answerdir / "rwi.nml"
    outdir = TEST_ARTIFACTS_DIR / stem
    itf = main(conffile, output_dir=outdir)
    answerfile = answerdir / "answer.pickle"
    # regold(itf, answerfile)


class TestRegressionVtuRef(AbstractTestImage):
    stem = "ref"
    answerdir = TEST_VTU_DATA_DIR / stem
    conffile = answerdir / "vac2fost_conf_nonaxisym.nml"
    outdir = TEST_ARTIFACTS_DIR /  f"vtu/{stem}"
    itf = main(conffile, output_dir=outdir, force_preroll=True, override={"flags": dict(read_gas_velocity=True)})
    answerfile = answerdir / "answer.pickle"
    # regold(itf, answerfile)


class TestRegressionDatRef(AbstractTestRegression):
    stem = "ref"
    answerdir = TEST_DAT_DATA_DIR / stem
    conffile = answerdir / "v2fconf.toml"
    outdir = TEST_ARTIFACTS_DIR / f"dat/{stem}"
    itf = main(conffile, output_dir=outdir, hydro_file_type="dat", force_preroll=True)
    answerfile = answerdir / "answer.pickle"
    # regold(itf, answerfile)


def test_dust_mass_estimations():
    itfs = [
        Interface(
            TEST_VTU_DATA_DIR / "vac2fost_conf_quick.nml",
            output_dir=TEST_ARTIFACTS_DIR,
            override={"flags": dict(dust_bin_mode=dbm)},
        )
        for dbm in ("mixed", "gas-only", "dust-only")
    ]
    for itf in itfs:
        itf.load_input_data()
    estimates = np.array([itf._estimate_dust_mass() for itf in itfs])
    np.testing.assert_allclose(estimates/estimates.min(), np.ones_like(estimates), rtol=1e-11)
    np.testing.assert_allclose(estimates, 4e-4, rtol=1e-3)