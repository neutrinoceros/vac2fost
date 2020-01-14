"""Commmon definitions for test files"""
from pathlib import Path

TEST_DIR = (Path(__file__).parent / "tests").resolve()
TEST_ARTIFACTS_DIR = TEST_DIR / ".artifacts"  # replace "output/"
TEST_DATA_DIR = TEST_DIR / "data"   # replace "sample/"
TEST_ANSWER_DIR = TEST_DATA_DIR / "answer"  # replace "ref/"
