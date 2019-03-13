from os import mkdir
from pathlib import Path
import pytest
from vac2fost.vac2fost import MCFOSTUtils
from vac2fost import main as app

testdir = Path(__file__).parent.parent
outdir = testdir/"output"

if not outdir.is_dir():
    mkdir(outdir)

def test_unicity():
    found = []
    for block, lines in MCFOSTUtils.blocks_descriptors.items():
        for line in lines:
            print("from test : ", block, line)
            found += [param for param in line]

    for p in found:
        if found.count(p) > 1:
            print(f"non unique key detected {p}")
    assert len(set(found)) == len(found)

def test_writter_null():
    MCFOSTUtils.write_mcfost_conf(outdir/"writter_out.para")

def test_writter_args():
    MCFOSTUtils.write_mcfost_conf(
        outdir/"writter_out_2.para",
        custom={"nphot_sed": 2}
    )


def test_unrecognized_mcfost_parameter():
    with pytest.raises(KeyError):
        app(
            str(testdir / "sample/vac2fost_conf_fake_params.nml"),
            output_dir=testdir/"output/fake_params"
        )
