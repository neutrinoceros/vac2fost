"""Commmon definitions for test files"""
from pathlib import Path
from shutil import rmtree

TEST_DIR = (Path(__file__).parent).resolve()
TEST_ARTIFACTS_DIR = TEST_DIR / ".artifacts"
TEST_DATA_DIR = TEST_DIR / "data"
TEST_VTU_DATA_DIR = TEST_DATA_DIR / "vtu"
TEST_DAT_DATA_DIR = TEST_DATA_DIR / "dat"

if TEST_ARTIFACTS_DIR.is_dir():
    rmtree(TEST_ARTIFACTS_DIR)