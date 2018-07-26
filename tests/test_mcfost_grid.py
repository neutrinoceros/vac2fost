'''Test basic consistency of main() called through python.'''

from pathlib import Path
import f90nml

from amrvac_pywrap import read_amrvac_conf
from vac2fost.vac2fost import MCFOSTUtils

here = Path(__file__).absolute().parent


config = f90nml.read(here / 'sample/vac2fost_conf.nml')
options = config['target_options']
sim_conf = read_amrvac_conf(files=options['amrvac_conf'], origin=options['origin'])

custom = {}
custom.update(MCFOSTUtils.translate_amrvac_conf(sim_conf))

def test_get_grid():
    output_dir = here / 'output/test_get_grid/'

    custom.update(config['mcfost_list'])
    mcfost_para_file = str(output_dir/'mcfost_conf.para')
    MCFOSTUtils.write_mcfost_conf(
        output_file=mcfost_para_file,
        custom=custom,
    )

    success = False
    try:
        target_grid = MCFOSTUtils.get_mcfost_grid(
            mcfost_conf=mcfost_para_file,
            mcfost_list=config['mcfost_list'],
            output_dir=output_dir,
            silent=False
        )
        success = True
    finally:
        assert success

def test_get_large_grid():
    output_dir = here / 'output/test_get_large_grid/'

    custom.update(config['mcfost_list'])
    mcfost_para_file = str(output_dir/'mcfost_conf.para')
    MCFOSTUtils.write_mcfost_conf(
        output_file=mcfost_para_file,
        custom=custom,
    )

    #overwrite
    config['mcfost_list']['nr'] = 200
    config['mcfost_list']['nphi'] = 200

    success = False
    try:
        target_grid = MCFOSTUtils.get_mcfost_grid(
            mcfost_conf=mcfost_para_file,
            mcfost_list=config['mcfost_list'],
            output_dir=output_dir,
            silent=False
        )
        success = True
    finally:
        assert success
