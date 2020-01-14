from os import mkdir
from shutil import rmtree
import multiprocessing as mp
import pytest

from vac2fost import main
from vac2fost.mcfost_utils import blocks_descriptors, write_mcfost_conf, get_mcfost_grid, get_mcfost_grid_dict
from vac2fost.logger import v2flogger as log

from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR


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
    write_mcfost_conf(TEST_ARTIFACTS_DIR / "writter_out.para")


def test_writter_args():
    write_mcfost_conf(TEST_ARTIFACTS_DIR / "writter_out_2.para", custom_parameters={"nphot_sed": 2})


def test_unrecognized_mcfost_parameter():
    with pytest.raises(ValueError):
        main(
            str(TEST_DATA_DIR / "vac2fost_conf_fake_params.nml"),
            output_dir=TEST_ARTIFACTS_DIR / "fake_params",
        )


def test_get_grid():
    output_dir = TEST_ARTIFACTS_DIR / "test_mcfost_api/"
    kwargs = dict(mcfost_conf_file=TEST_DATA_DIR/"ref3.0.para",
                output_dir=output_dir / "dry_grid", require_run=True)
    grid = get_mcfost_grid(**kwargs)
    grid_dict = get_mcfost_grid_dict(**kwargs)
