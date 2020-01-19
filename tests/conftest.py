"""Commmon definitions for test files"""
from pathlib import Path
from shutil import rmtree

TEST_DIR = (Path(__file__).parent).resolve()
TEST_ARTIFACTS_DIR = TEST_DIR / ".artifacts"  # replace "output/"
TEST_DATA_DIR = TEST_DIR / "data"   # replace "sample/"
TEST_ANSWER_DIR = TEST_DATA_DIR / "answer"  # replace "ref/"

if TEST_ARTIFACTS_DIR.is_dir():
    rmtree(TEST_ARTIFACTS_DIR)