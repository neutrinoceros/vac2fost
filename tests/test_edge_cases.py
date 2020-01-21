from os import chdir
from shutil import rmtree
from subprocess import run
from astropy.io import fits
import f90nml

from vac2fost import main
from vac2fost.logger import v2flogger
from vac2fost.interfaces import VtuFileInterface

import pytest
from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR
densfile = "hd142527_dusty0000.fits"


def test_init_from_toml():
    try:
        import toml
        VtuFileInterface(TEST_DATA_DIR / "vac2fost_conf_base.toml")
    except ImportError:
        pass


@pytest.mark.parametrize("mode,expected_n_hdu",
                        [("gas-only", 2), ("dust-only", 3), ("mixed", 3)])
def test_read_gas_vel_with_dbm_mode(mode, expected_n_hdu):
    outdir = TEST_ARTIFACTS_DIR / f"gasvel_{mode}"
    if outdir.exists():
        rmtree(outdir)
    main(
        TEST_DATA_DIR / "vac2fost_conf_quick.nml",
        override={"flags": dict(read_gas_velocity=True, dust_bin_mode=mode)},
        output_dir=outdir,
    )
    chdir(outdir)
    with fits.open(densfile) as dfile:
        assert len(dfile) == expected_n_hdu
        dat = dfile[-1].data
        assert dat.shape == (3, 10, 2, 10)

    # check mcfost can still be run
    run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)


# @pytest.mark.incremental #each test is run only if the previous one passed
class TestReadGasDensity:
    outdir = TEST_ARTIFACTS_DIR / "test_read_gas_density"
    if outdir.is_dir():
        rmtree(outdir)
    conf_file = TEST_DATA_DIR / "vac2fost_conf_quick.nml"
    v2flogger.setLevel(10)
    itf = main(
        conf_file,
        output_dir=outdir,
        override={"flags": dict(dust_bin_mode="dust-only", read_gas_density=True)},
    )

    def test_read_gas_density_header(self):
        header = fits.open(__class__.itf.io.OUT.filepath)[0].header
        cards_as_dicts = {c.keyword: c.value for c in header.cards}
        assert ("read_gas_density", 1) in cards_as_dicts.items()

    def test_read_gas_density_shape(self):
        gas_density = fits.open(__class__.itf.io.OUT.filepath)[-1].data
        conf = f90nml.read(__class__.conf_file)["mcfost_output"]
        target_shape = tuple([conf[k] for k in ("n_rad", "n_z", "n_az")])
        assert gas_density.shape == target_shape

    def test_read_gas_density_is_valid(self):
        """check that mcfost doesn't crash when passed gas density."""
        chdir(__class__.itf.io.OUT.directory)
        run(
            f"mcfost mcfost_conf.para -density_file {__class__.itf.io.OUT.filepath.name} -3D",
            shell=True,
            check=True,
        )
