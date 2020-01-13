import pathlib
import numpy as np
import pytest
import f90nml
from vtk_vacreader import VacDataSorter as VDS
from vac2fost import VtuFileInterface as Interface
from vac2fost.mcfost_utils import MINGRAINSIZE_mum

test_dir = pathlib.Path(__file__).absolute().parent.parent
output_dir = test_dir / "output/test_dbm"


class TestDBM:
    """test automatic behaviour of dust binning mode selection"""

    def test_unrecognized_dbm(self):
        with pytest.raises(ValueError):
            Interface(
                test_dir / "sample/vac2fost_conf_quick_no_dust.nml",
                output_dir=output_dir,
                override={"flags": dict(dust_bin_mode="not-a-real-dbm-option")},
            )

    @pytest.mark.parametrize("mode,expected",
        [("gas-only", [MINGRAINSIZE_mum]),
         ("dust-only", [MINGRAINSIZE_mum, 1e4, 1e3]),
         ("mixed", [MINGRAINSIZE_mum, 1e4, 1e3])])
    def test_grain_sizes_with_mode(self, mode, expected):
        itf = Interface(
            test_dir / "sample/vac2fost_conf_quick.nml",
            output_dir=output_dir,
            override={"flags": dict(dust_bin_mode=mode)},
        )
        assert itf._dust_bin_mode == mode
        np.testing.assert_array_equal(itf._grain_micron_sizes, expected)

    def test_auto_into_mixed(self):
        itf = Interface(test_dir / "sample/vac2fost_conf_quick.nml", output_dir=output_dir)
        assert itf._dust_bin_mode == "mixed"

    def test_auto_into_gas_only(self):
        itf = Interface(test_dir / "sample/vac2fost_conf_quick_no_dust.nml", output_dir=output_dir)
        assert itf._dust_bin_mode == "gas-only"

    @pytest.mark.parametrize("mode", ["mixed", "dust-only"])
    def test_impossible_mode_into_failure(self, mode):
        with pytest.raises(KeyError):
            Interface(
                test_dir / "sample/vac2fost_conf_quick_no_dust.nml",
                output_dir=output_dir,
                override={"flags": dict(dust_bin_mode=mode)},
            )


class TestMassEstimate:
    def test_dust_mass_estimations_consistency(self):
        itfs = [
            Interface(
                test_dir / "sample/vac2fost_conf_quick.nml",
                output_dir=output_dir,
                override={"flags": dict(dust_bin_mode=dbm)},
            )
            for dbm in ("mixed", "gas-only", "dust-only")
        ]
        for itf in itfs:
            itf.load_input_data()
        estimates = [itf._estimate_dust_mass() for itf in itfs]
        ref = estimates.pop(0)
        for e in estimates:
            assert abs(e - ref) / ref < 1e-11

    def test_dust_mass_absolute_estimation(self):
        VDS(str(test_dir / "sample/flat_rphi0000.vtu"), shape=(512, 128))
        conf = f90nml.read(test_dir / "sample/flat_rphi.nml")
        rmin = conf["meshlist"]["xprobmin1"]
        rmax = conf["meshlist"]["xprobmax1"]
        sig0 = conf["disk_list"]["rho0"]
        g2d = conf["usr_dust_list"]["gas2dust_ratio"]
        ref = np.pi * sig0 * (rmax ** 2 - rmin ** 2) / g2d
        for dbm in ("mixed", "gas-only", "dust-only"):
            itf = Interface(
                test_dir / "sample/vac2fost_conf_flatdisk.nml",
                override={"flags": dict(dust_bin_mode=dbm)},
            )
            itf.load_input_data()
            assert abs(itf._estimate_dust_mass() - ref) / ref < 1e-7
