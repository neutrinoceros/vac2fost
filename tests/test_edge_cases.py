from os import chdir
from subprocess import run
from pathlib import Path
from vac2fost import main as app

testdir = Path(__file__).parent.resolve()
OUT = testdir/"output"

def test_read_gas_vel_gas_only():
    outdir = OUT / "gasvel_gasonly"
    app(
        testdir / "sample/vac2fost_conf_quick.nml",
        read_gas_velocity=True,
        dust_bin_mode="gas-only",
        output_dir=outdir
    )
    chdir(outdir)
    run(["mcfost", "mcfost_conf.para", "-3D", "-density_file", "hd142527_dusty0000.fits"], check=True)
