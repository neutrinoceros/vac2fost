#!/usr/bin/env python3
"""
A convenient reusable script embedding most of the package's functionality.

For documentation on usage, run `python v2f.py --help`

The main algorithm runs the following steps
  a) load AMRVAC data with vtk_vacreader.VacDataSorter(), sort it as 2D arrays
  b) dry-run MCFOST to get the exact output
  c) interpolate data to the target grid
  d) convert to 3D (gaussian redistribution of density)
  e) collect, sort and write output data to a fits file


Known limitations
  1) AMR grids are not supported (.vtu files need to be converted to uniform grids)
  2) .vtu are assumed to use polar coordinates (r-phi 2D)
  3) 2D interpolation does not account for the curvature of polar cells
"""
import sys
from pathlib import Path
from argparse import ArgumentParser

import f90nml

from vac2fost.info import __version__
from vac2fost.mcfost_utils import blocks_descriptors
from vac2fost.utils import CYAN, BOLD, get_prompt_size, decorated_centered_message
from vac2fost.interfaces import DEFAULT_UNITS, Interface, VerbatimInterface
from vac2fost.main import main


# Definitions =======================================================================
def generate_conf_template() -> f90nml.Namelist:
    """Generate a template namelist object with comments instead of default values"""
    amrvac_list = dict(
        hydro_data_dir="path/to/output/data/directory",
        config="relative/to/<hydro_data_dir>/path/to/amrvac/config/file[s]",
        nums=0
    )

    mcfost_list = dict(
        n_rad=128, n_rad_in=4, n_az=128, nz=10,
        # aspect ratio is implied by those parameters
        flaring_exp=1.125,
        reference_radius=100.0,  # [a.u.]
        scale_height=1.0,  # [a.u.], at defined at reference_radius
    )
    sublists = {
        "amrvac_input": amrvac_list,
        "units": DEFAULT_UNITS,
        "mcfost_output": mcfost_list
    }
    template = f90nml.Namelist({k: f90nml.Namelist(v) for k, v in sublists.items()})
    return template


def print_mcfost_default_conf():
    from pprint import pprint
    for block_name, lines in blocks_descriptors.items():
        print(block_name)
        print("="*len(block_name))
        for line in lines:
            for k, v in line.items():
                print(f"{k} = {v}")
        print()

# Script ============================================================================
if __name__ == "__main__":
    # Parse the script arguments
    parser = ArgumentParser(description="Parse arguments for main app")
    parser.add_argument(
        dest="configuration", type=str,
        nargs="?",
        default=None,
        help="configuration file (namelist) for this script"
    )
    parser.add_argument(
        "-n", "--nums", dest="nums", type=int,
        required=False,
        default=None,
        nargs="*",
        help="output number(s) of the target .vtu VAC output file to be converted"
    )
    parser.add_argument(
        "-o", "--output", dest="output", type=str,
        required=False,
        default=".",
        help="select output directory for generated files"
    )
    parser.add_argument(
        "-dbm", "--dustbinmode", dest="dbm", type=str,
        required=False,
        default="auto",
        help="prefered bin selection mode [dust-only, gas-only, mixed, auto]"
    )
    parser.add_argument(
        "--read_gas_density",
        action="store_true",
        help="pass gas density to mcfost"
    )
    parser.add_argument(
        "--settling",
        action="store_true",
        help="Dubrulle-style grain settling"
    )
    parser.add_argument(
        "--read_gas_velocity",
        action="store_true",
        help="pass gas velocity to mcfost (keplerian velocity is assumed otherwise)"
    )
    parser.add_argument(
        "--axisymmetry",
        action="store_true",
        help="generate a r-z slice"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="verbose logging"
    )
    parser.add_argument(
        "--mcfost_verbose",
        action="store_true",
        help="do not silence mcfost"
    )
    parser.add_argument(
        "--genconf", action="store_true",
        help="print a default configuration file for vac2fost"
    )
    parser.add_argument(
        "--print_mcfost_defaults", action="store_true",
        help="print internal default values for all mcfost parameters as returned by vac2fost"
    )
    parser.add_argument(
        "--cprofile",
        action="store_true",
        help="activate code profiling"
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="display vac2fost version"
    )

    cargs = parser.parse_args()

    if cargs.version:
        print(__version__)
        sys.exit(0)

    # special cases -----------------------------
    if cargs.genconf:
        print(generate_conf_template())
        print(f"%% automatically generated with vac2fost {__version__}\n")
        sys.exit(0)
    elif cargs.print_mcfost_defaults:
        print_mcfost_default_conf()
        sys.exit(0)
    elif len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if cargs.cprofile:
        import cProfile
        import pstats
        import io
        pr = cProfile.Profile()
        pr.enable()
    # -------------------------------------------
    main(
        config_file=cargs.configuration,
        nums=cargs.nums,
        output_dir=cargs.output.strip(),
        dust_bin_mode=cargs.dbm,
        read_gas_density=cargs.read_gas_density,
        read_gas_velocity=cargs.read_gas_velocity,
        settling=cargs.settling,
        axisymmetry=cargs.axisymmetry,
        loglevel={True: 0, False: 30}[cargs.verbose],
        mcfost_verbose=cargs.mcfost_verbose # wip
    )
    # -------------------------------------------
    if cargs.cprofile:
        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
        ps.print_stats()
        print(s.getvalue())
