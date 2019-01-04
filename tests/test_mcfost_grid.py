'''Test basic consistency of main() called through python.'''
import shutil
import multiprocessing as mp
from os import mkdir
from pathlib import Path

import pytest
import f90nml

from amrvac_pywrap import read_amrvac_conf
from vac2fost.vac2fost import MCFOSTUtils, Interface

test_dir = Path(__file__).absolute().parent

sampledir = test_dir / 'sample'
config_file = sampledir / 'vac2fost_conf.nml'
config = f90nml.read(config_file)
options = config['target_options']
sim_conf = read_amrvac_conf(files=options['amrvac_conf'], origin=sampledir/options['origin'])

custom = {}
itf = Interface(str(config_file))
custom.update(MCFOSTUtils.translate_amrvac_conf(itf))

def gen_mcfost_grid(output_dir):
    """Create an Interface only to create an mcfost grid from it,
    this is intented to be called in parallel
    """
    if output_dir.exists():
        shutil.rmtree(output_dir)
    itf = Interface(config_file=test_dir/'sample/vac2fost_conf.nml', output_dir=output_dir, dbg=True)
    return itf.output_grid

def test_dry_grid_gen():
    """Check that the output grid can be retrieved simply by calling it at the interface level."""
    output_dir = test_dir/"output/test_parallel_instanciation/"
    gen_mcfost_grid(output_dir/"dry_grid")

@pytest.mark.skipif(mp.cpu_count()==1, reason="parallel computation only with Ncpus>=2")
def test_parallel_instanciation():
    """Check that instanciating multiple Interface class object at same
    time and location doesn't create collisions"""
    output_dir = test_dir/"output/test_parallel_instanciation/"
    ncpu = min(mp.cpu_count(), 4)
    with mp.Pool(ncpu) as pool:
        pool.map(gen_mcfost_grid, [output_dir/str(i) for i in range(ncpu)])

def test_get_grid():
    output_dir = test_dir / 'output/test_get_grid/'
    if not output_dir.is_dir():
        mkdir(output_dir)
    custom.update(config['mcfost_list'])
    mcfost_para_file = str(output_dir/'mcfost_conf.para')
    MCFOSTUtils.write_mcfost_conf(
        output_file=mcfost_para_file,
        custom=custom,
    )

    success = False
    try:
        target_grid = MCFOSTUtils.get_mcfost_grid(
            mcfost_conf_file=mcfost_para_file,
            mcfost_list=config['mcfost_list'],
            output_dir=output_dir,
            silent=False
        )
        success = True
    finally:
        assert success

def test_get_large_grid():
    output_dir = test_dir / 'output/test_get_large_grid/'
    if not output_dir.is_dir():
        mkdir(output_dir)

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
            mcfost_conf_file=mcfost_para_file,
            mcfost_list=config['mcfost_list'],
            output_dir=output_dir,
            silent=False
        )
        success = True
    finally:
        assert success
