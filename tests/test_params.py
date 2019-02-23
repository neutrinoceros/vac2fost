from pathlib import Path
import os

from numpy import pi
import pytest
import f90nml

from vtk_vacreader import VacDataSorter as VDS
from vac2fost.vac2fost import get_dust_mass
from vac2fost import main as app

here = Path(__file__).absolute().parent
outdir = here / 'output/fake_params'
gridfile = outdir / 'mcfost_grid.fits.gz'
if gridfile.exists():
    os.remove(gridfile)

def test_dust_mass_estimation():
    data = VDS(str(here/'sample/flat_rphi0000.vtu'), shape=(512, 128))
    estimate = get_dust_mass(data)
    conf = f90nml.read(here / 'sample/flat_rphi.nml')
    rmin = conf['meshlist']['xprobmin1']
    rmax = conf['meshlist']['xprobmax1']
    sig0 = conf['disk_list']['rho0']
    g2d = conf['usr_dust_list']['gas2dust_ratio']
    ref = pi * sig0 * (rmax**2 - rmin**2) / g2d
    err = abs(estimate - ref) / ref
    assert err < 1e-7

def test_unrecognized_mcfost_parameter():
    with pytest.raises(KeyError):
        app(
            str(here/'sample/vac2fost_conf_fake_params.nml'),
            output_dir=outdir
        )
