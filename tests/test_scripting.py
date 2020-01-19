"""Test basic consistency of main() called through python."""
import shutil

import f90nml

from vac2fost import main
from vac2fost.logger import v2flogger as log

import pytest
from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR

log.setLevel(10)
output_dir = TEST_ARTIFACTS_DIR / "test_py_scripting"
if output_dir.exists():
    shutil.rmtree(output_dir)

def test_python_call():
    main(TEST_DATA_DIR / "vac2fost_conf.nml", output_dir=output_dir)

def test_python_call_multiple():
    """Note: we purposedly ask to convert data that is not available (n=3 and n=9999)
    to check that the main function handles it correctly (warning + resume iteration)"""
    nums = [1, 3, 9999]
    itf = main(
        TEST_DATA_DIR / "vac2fost_conf.nml",
        output_dir=output_dir,
        override={"amrvac_input": {"nums": nums}},
        loglevel=10
    )
    assert (output_dir / f"hd142527_dusty0001.fits").exists()

    # check that the internal iterator was fully consumed
    with pytest.raises(StopIteration):
        itf.advance_iteration()
