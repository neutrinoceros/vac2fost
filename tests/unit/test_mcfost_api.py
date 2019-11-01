from os import mkdir
from shutil import rmtree
from pathlib import Path
import multiprocessing as mp
import pytest

from vac2fost.mcfost_utils import blocks_descriptors, write_mcfost_conf
from vac2fost import Interface, main as app

testdir = Path(__file__).parent.parent
outdir = testdir/"output"

if not outdir.is_dir():
    mkdir(outdir)

def test_unicity():
    found = []
    for block, lines in blocks_descriptors.items():
        for line in lines:
            print("from test : ", block, line)
            found += [param for param in line]

    for p in found:
        if found.count(p) > 1:
            print(f"non unique key detected {p}")
    assert len(set(found)) == len(found)

def test_writter_null():
    write_mcfost_conf(outdir/"writter_out.para")

def test_writter_args():
    write_mcfost_conf(
        outdir/"writter_out_2.para",
        custom_parameters={"nphot_sed": 2}
    )


def test_unrecognized_mcfost_parameter():
    with pytest.raises(KeyError):
        app(
            str(testdir / "sample/vac2fost_conf_fake_params.nml"),
            output_dir=testdir/"output/fake_params"
        )


def gen_mcfost_grid(output_dir):
    """Create an Interface only to create an mcfost grid from it,
    this is intented to be called in parallel
    """
    if output_dir.exists():
        rmtree(output_dir)
    myitf = Interface(config_file=testdir / "sample/vac2fost_conf.nml",
                      output_dir=output_dir, mcfost_verbose=True)
    return myitf.output_grid

class TestMcfostGridGen:
    def test_dry_grid_gen(self):
        """Check that the output grid can be retrieved simply by calling it at the interface level."""
        output_dir = testdir/"output/test_parallel_instanciation/"
        if output_dir.exists():
            rmtree(output_dir)
        gen_mcfost_grid(output_dir/"dry_grid")

    @pytest.mark.skipif(mp.cpu_count() == 1, reason="parallel computation only with Ncpus>=2")
    def test_parallel_instanciation(self):
        """Check that instanciating multiple Interface class object at same
        time and location doesn't create collisions"""
        output_dir = testdir/"output/test_parallel_instanciation/"
        if output_dir.exists():
            rmtree(output_dir)
        ncpu = min(mp.cpu_count(), 4)
        with mp.Pool(ncpu) as pool:
            pool.map(gen_mcfost_grid, [output_dir/str(i) for i in range(ncpu)])
