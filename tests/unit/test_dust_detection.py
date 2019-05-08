import pathlib
import numpy as np
import pytest
import f90nml
from vtk_vacreader import VacDataSorter as VDS
from vac2fost import Interface
from vac2fost.mcfost_utils import MINGRAINSIZE_µ

test_dir = pathlib.Path(__file__).absolute().parent.parent
output_dir = test_dir / "output/test_dbm"

class TestDBM:
    """test automatic behaviour of dust binning mode selection"""
    def test_unrecognized_dbm(self):
        with pytest.raises(KeyError):
            Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                      output_dir=output_dir,
                      dust_bin_mode="not-a-real-dbm-option")

    def test_gas_only(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                        output_dir=output_dir,
                        dust_bin_mode="gas-only")
        assert itf.dust_binning_mode == "gas-only"
        assert itf.grain_micron_sizes == [MINGRAINSIZE_µ]

    def test_dust_only(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                        output_dir=output_dir,
                        dust_bin_mode="dust-only")
        assert itf.dust_binning_mode == "dust-only"
        np.testing.assert_array_equal(itf.grain_micron_sizes, [MINGRAINSIZE_µ, 1e4, 1e3])

    def test_mixed(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                        output_dir=output_dir,
                        dust_bin_mode="mixed")
        assert itf.dust_binning_mode == "mixed"
        np.testing.assert_array_equal(itf.grain_micron_sizes, [MINGRAINSIZE_µ, 1e4, 1e3])

    def test_auto_into_mixed(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                        output_dir=output_dir,
                        dust_bin_mode="auto")
        assert itf.dust_binning_mode == "mixed"

    def test_auto_into_gas_only(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                        output_dir=output_dir,
                        dust_bin_mode="auto")
        assert itf.dust_binning_mode == "gas-only"

    def test_mixed_into_KeyError(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                        output_dir=output_dir,
                        dust_bin_mode="mixed")
        with pytest.raises(KeyError):
            # impossible setups should raise KeyError
            gms = itf.grain_micron_sizes

    def test_dust_only_into_KeyError(self):
        itf = Interface(test_dir/'sample/vac2fost_conf_quick_no_dust.nml',
                        output_dir=output_dir,
                        dust_bin_mode="dust-only")
        with pytest.raises(KeyError):
            gms = itf.grain_micron_sizes



class TestMassEstimate:
    def test_dust_mass_estimations_consistency(self):
        itfs = [Interface(test_dir/'sample/vac2fost_conf_quick.nml',
                          output_dir=output_dir,
                          dust_bin_mode=dbm)
                for dbm in ("mixed", "gas-only", "dust-only")]
        estimates = [itf.estimate_dust_mass() for itf in itfs]
        ref = estimates.pop(0)
        for e in estimates:
            assert abs(e -ref) / ref < 1e-11

    def test_dust_mass_absolute_estimation(self):
        data = VDS(str(test_dir / "sample/flat_rphi0000.vtu"), shape=(512, 128))
        conf = f90nml.read(test_dir / "sample/flat_rphi.nml")
        rmin = conf["meshlist"]["xprobmin1"]
        rmax = conf["meshlist"]["xprobmax1"]
        sig0 = conf["disk_list"]["rho0"]
        g2d = conf["usr_dust_list"]["gas2dust_ratio"]
        ref = np.pi * sig0 * (rmax**2 - rmin**2) / g2d
        for dbm in ("mixed", "gas-only", "dust-only"):
            itf = Interface(test_dir/"sample/vac2fost_conf_flatdisk.nml", dust_bin_mode=dbm)
            assert abs(itf.estimate_dust_mass() - ref) / ref < 1e-7
