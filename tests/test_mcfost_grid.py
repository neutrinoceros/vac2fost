import shutil
import multiprocessing as mp
from pathlib import Path

import pytest

from vac2fost.vac2fost import Interface

test_dir = Path(__file__).absolute().parent

def gen_mcfost_grid(output_dir):
    """Create an Interface only to create an mcfost grid from it,
    this is intented to be called in parallel
    """
    if output_dir.exists():
        shutil.rmtree(output_dir)
    myitf = Interface(config_file=test_dir / "sample/vac2fost_conf.nml",
                      output_dir=output_dir, mcfost_verbose=True)
    return myitf.output_grid

def test_dry_grid_gen():
    """Check that the output grid can be retrieved simply by calling it at the interface level."""
    output_dir = test_dir/"output/test_parallel_instanciation/"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    gen_mcfost_grid(output_dir/"dry_grid")

@pytest.mark.skipif(mp.cpu_count() == 1, reason="parallel computation only with Ncpus>=2")
def test_parallel_instanciation():
    """Check that instanciating multiple Interface class object at same
    time and location doesn't create collisions"""
    output_dir = test_dir/"output/test_parallel_instanciation/"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    ncpu = min(mp.cpu_count(), 4)
    with mp.Pool(ncpu) as pool:
        pool.map(gen_mcfost_grid, [output_dir/str(i) for i in range(ncpu)])
