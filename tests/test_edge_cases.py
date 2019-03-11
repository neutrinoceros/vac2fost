from os import chdir
from shutil import rmtree
from subprocess import run
from pathlib import Path
import pytest
from astropy.io import fits
from vac2fost import main as app
from vac2fost.vac2fost import DETECTED_MCFOST_VERSION

testdir = Path(__file__).parent.resolve()
OUT = testdir/"output"
densfile = "hd142527_dusty0000.fits"


@pytest.mark.skipif(DETECTED_MCFOST_VERSION < "3.0.35", reason="bug only with mcfost 3.0.35")
def test_read_gas_vel_gas_only():
    outdir = OUT / "gasvel_gasonly"
    if outdir.exists():
        rmtree(outdir)
    app(
        testdir / "sample/vac2fost_conf_quick.nml",
        read_gas_velocity=True,
        dust_bin_mode="gas-only",
        output_dir=outdir
    )
    chdir(outdir)
    with fits.open(densfile) as dfile:
        dat = dfile[2].data
        assert dat.shape == (3, 10, 2, 10)

    run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)

@pytest.mark.skipif(DETECTED_MCFOST_VERSION < "3.0.35", reason="bug only with mcfost 3.0.35")
def test_read_gas_vel_dust_only():
    outdir = OUT / "gasvel_dustonly"
    if outdir.exists():
        rmtree(outdir)
    app(
        testdir / "sample/vac2fost_conf_quick.nml",
        read_gas_velocity=True,
        dust_bin_mode="dust-only",
        output_dir=outdir
    )
    chdir(outdir)
    with fits.open(densfile) as dfile:
        dat = dfile[2].data
        assert dat.shape == (3, 10, 2, 10)

    run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)

@pytest.mark.skipif(DETECTED_MCFOST_VERSION < "3.0.35", reason="bug only with mcfost 3.0.35")
def test_read_gas_vel_mixed():
    outdir = OUT / "gasvel_mixed"
    if outdir.exists():
        rmtree(outdir)
    app(
        testdir / "sample/vac2fost_conf_quick.nml",
        read_gas_velocity=True,
        dust_bin_mode="mixed",
        output_dir=outdir
    )
    chdir(outdir)
    with fits.open(densfile) as dfile:
        dat = dfile[2].data
        assert dat.shape == (3, 10, 2, 10)

    run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", densfile], check=True)
