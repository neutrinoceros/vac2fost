'''Test basic consistency of main() called through python.'''

from pathlib import Path
import f90nml

from vac2fost.vac2fost import MCFOSTUtils

here = Path(__file__).absolute().parent

def test_get_grid():
    config = f90nml.read(here / 'sample/vac2fost_conf.nml')

    mesh_dimensions = {'xprobmin1': 70, 'xprobmax1': 450}
    success = False
    try:
        target_grid = MCFOSTUtils.get_mcfost_grid(
            mcfost_list=config['mcfost_list'],
            mesh_list=mesh_dimensions,
            silent=False
        )
        success = True
    finally:
        assert success

def test_get_large_grid():
    config = f90nml.read(here / 'sample/vac2fost_conf.nml')

    config['mcfost_list']['nr'] = 1024
    config['mcfost_list']['nphi'] = 1024

    mesh_dimensions = {'xprobmin1': 70, 'xprobmax1': 450}
    success = False
    try:
        target_grid = MCFOSTUtils.get_mcfost_grid(
            mcfost_list=config['mcfost_list'],
            mesh_list=mesh_dimensions,
            silent=False
        )
        success = True
    finally:
        assert success
