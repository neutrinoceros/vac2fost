import numpy as np
import pytest
import f90nml
from vtk_vacreader import VacDataSorter as VDS
from vac2fost import VtuFileInterface as Interface
from vac2fost.mcfost_utils import MINGRAINSIZE_mum

from conftest import TEST_DATA_DIR, TEST_ARTIFACTS_DIR

output_dir = TEST_ARTIFACTS_DIR / "test_dust_detection"

class TestDBM:
    """test automatic behaviour of dust binning mode selection"""

    def test_unrecognized_dbm(self):
        with pytest.raises(ValueError):
            Interface(
                TEST_DATA_DIR / "vac2fost_conf_quick_no_dust.nml",
                output_dir=output_dir,
                override={"flags": dict(dust_bin_mode="not-a-real-dbm-option")},
            )

    @pytest.mark.parametrize("mode,expected",
                             [("gas-only", [MINGRAINSIZE_mum]),
                             ("dust-only", [MINGRAINSIZE_mum, 1e4, 1e3]),
                             ("mixed", [MINGRAINSIZE_mum, 1e4, 1e3])])
    def test_grain_sizes_with_mode(self, mode, expected):
        itf = Interface(
            TEST_DATA_DIR / "vac2fost_conf_quick.nml",
            output_dir=output_dir,
            override={"flags": dict(dust_bin_mode=mode)},
        )
        assert itf._dust_bin_mode == mode
        np.testing.assert_array_equal(itf._grain_micron_sizes, expected)

    @pytest.mark.parametrize("conf_file,expected",
                             [("vac2fost_conf_quick.nml", "mixed"),
                              ("vac2fost_conf_quick_no_dust.nml", "gas-only")])
    def test_auto_dbm(self, conf_file, expected):
        itf = Interface(TEST_DATA_DIR / conf_file, output_dir=output_dir)
        assert itf._dust_bin_mode == expected

    @pytest.mark.parametrize("mode", ["mixed", "dust-only"])
    def test_impossible_mode_into_failure(self, mode):
        with pytest.raises(KeyError):
            Interface(
                TEST_DATA_DIR / "vac2fost_conf_quick_no_dust.nml",
                output_dir=output_dir,
                override={"flags": dict(dust_bin_mode=mode)},
            )
