"""Tests related to isolated functions in vac2fost (outside of classes)"""
import os
from pathlib import Path
from vac2fost.utils import shell_path


def test_expandvars():
    os.environ["VAC2FOST_FAKE"] = "/fake/path/to/nothing"
    assert shell_path("$VAC2FOST_FAKE") == Path(os.path.expandvars("$VAC2FOST_FAKE"))

def test_expanduser():
    assert shell_path("~") == Path("~").expanduser()
    assert shell_path("$HOME") == Path.home()
