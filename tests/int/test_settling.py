"""Test on grain-size differentiate settling mode"""
from shutil import rmtree
from pathlib import Path

from vac2fost.logger import v2flogger as log
from vac2fost import main as app


testdir = Path(__file__).parent.parent
OUT = testdir / "output"
densfile = "hd142527_dusty0000.fits"


class TestSettling:
    outdir = testdir / "output/test_settling"
    if outdir.is_dir():
        rmtree(outdir)
    conf_file = testdir / "sample/vac2fost_conf_quick.nml"

    def test_default_settling(self):
        log.setLevel(10)
        itf = app(__class__.conf_file, output_dir=__class__.outdir)
        assert not itf._use_settling

    def test_activate_settling(self):
        log.setLevel(10)
        itf = app(
            __class__.conf_file,
            output_dir=__class__.outdir,
            override={"flags": dict(settling=True)},
        )
        assert itf._use_settling

    def test_run_settling(self):
        log.setLevel(10)
        itf = app(
            __class__.conf_file,
            output_dir=__class__.outdir,
            override={"flags": dict(settling=True)},
        )
        itf.preroll_mcfost()
        itf.get_output_ndarray()
