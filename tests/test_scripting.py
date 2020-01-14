"""Test basic consistency of main() called through python."""
import shutil

import f90nml

from vac2fost import main
from vac2fost.logger import v2flogger as log

from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR

log.setLevel(10)
output_dir = TEST_ARTIFACTS_DIR / "test_py_scripting"
if output_dir.exists():
    shutil.rmtree(output_dir)

def test_python_call():
    main(TEST_DATA_DIR / "vac2fost_conf.nml", output_dir=output_dir)

def test_python_call_multiple():
    nums = [0, 1]
    main(
        TEST_DATA_DIR / "vac2fost_conf.nml",
        output_dir=output_dir,
        override={"amrvac_input": {"nums": nums}},
    )
    for n in nums:
        assert (output_dir / f"hd142527_dusty{str(n).zfill(4)}.fits").exists()
