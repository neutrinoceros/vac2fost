from os import chdir
from shutil import rmtree
from subprocess import run
from pathlib import Path
from astropy.io import fits
import f90nml
from vac2fost import main as app
from vac2fost.logger import v2flogger

testdir = Path(__file__).parent.parent
OUT = testdir / "output"
densfile = "hd142527_dusty0000.fits"


class TestReadGasVelocity:
    def test_read_gas_vel_gas_only(self):
        outdir = OUT / "gasvel_gasonly"
        if outdir.exists():
            rmtree(outdir)
        app(
            testdir / "sample/vac2fost_conf_quick.nml",
            override={"flags": dict(read_gas_velocity=True, dust_bin_mode="gas-only")},
            output_dir=outdir,
        )
        chdir(outdir)
        with fits.open(densfile) as dfile:
            assert len(dfile) == 2
            dat = dfile[1].data
            assert dat.shape == (3, 10, 2, 10)

        run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)

    def test_read_gas_vel_dust_only(self):
        outdir = OUT / "gasvel_dustonly"
        if outdir.exists():
            rmtree(outdir)
        app(
            testdir / "sample/vac2fost_conf_quick.nml",
            override={"flags": dict(read_gas_velocity=True, dust_bin_mode="dust-only")},
            output_dir=outdir,
        )
        chdir(outdir)
        with fits.open(densfile) as dfile:
            assert len(dfile) == 3
            dat = dfile[2].data
            assert dat.shape == (3, 10, 2, 10)

        run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)

    def test_read_gas_vel_mixed(self):
        outdir = OUT / "gasvel_mixed"
        if outdir.exists():
            rmtree(outdir)
        app(
            testdir / "sample/vac2fost_conf_quick.nml",
            override={"flags": dict(read_gas_velocity=True, dust_bin_mode="mixed")},
            output_dir=outdir,
        )
        chdir(outdir)
        with fits.open(densfile) as dfile:
            assert len(dfile) == 3
            dat = dfile[2].data
            assert dat.shape == (3, 10, 2, 10)

        run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)


# @pytest.mark.incremental #each test is run only if the previous one passed
class TestReadGasDensity:
    outdir = testdir / "output/test_read_gas_density"
    if outdir.is_dir():
        rmtree(outdir)
    conf_file = testdir / "sample/vac2fost_conf_quick.nml"
    v2flogger.setLevel(10)
    itf = app(
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
        target_shape = tuple([conf[k] for k in ("n_rad", "nz", "n_az")])
        assert gas_density.shape == target_shape

    def test_read_gas_density_is_valid(self):
        """check that mcfost doesn't crash when passed gas density."""
        chdir(__class__.itf.io.OUT.directory)
        run(
            f"mcfost mcfost_conf.para -density_file {__class__.itf.io.OUT.filepath.name} -3D",
            shell=True,
            check=True,
        )
