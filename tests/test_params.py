from pathlib import Path
import os

from numpy import pi
import pytest
import f90nml

from vac2fost import main as app

here = Path(__file__).absolute().parent
outdir = here / 'output/fake_params'
gridfile = outdir / 'mcfost_grid.fits.gz'
if gridfile.exists():
    os.remove(gridfile)

def test_unrecognized_mcfost_parameter():
    with pytest.raises(KeyError):
        app(
            str(here/'sample/vac2fost_conf_fake_params.nml'),
            output_dir=outdir
        )
