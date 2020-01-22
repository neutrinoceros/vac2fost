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
from argparse import ArgumentParser

import f90nml

from vac2fost.info import __version__
from vac2fost.mcfost_utils import blocks_descriptors
from vac2fost.interfaces import DEFAULT_UNITS
from vac2fost.main import main


# Definitions =======================================================================
amrvac_list = dict(
    hydro_data_dir="path/to/output/data/directory",
    config="relative/to/<hydro_data_dir>/path/to/amrvac/config/file[s]",
    nums=0,
)
mcfost_list = dict(
    n_rad=128,
    n_rad_in=4,
    n_az=128,
    n_z=10,
    # aspect ratio is implied by those parameters
    flaring_exp=1.125,
    reference_radius=100.0,
    scale_height=1.0,  # [a.u.], at defined at reference_radius
)
dust_list = dict(
    grain_size2micron=1e4, grain_sizes=[1, 10, 100], dust_to_gas_ratio=0.01  # using original unit
)
sublists = {
    "amrvac_input": amrvac_list,
    "mcfost_output": mcfost_list,
    "dust": dust_list,
    "units": DEFAULT_UNITS,
    "flags": {},
}
template = f90nml.Namelist({k: f90nml.Namelist(v) for k, v in sublists.items()})


def print_mcfost_default_conf():
    for block_name, lines in blocks_descriptors.items():
        print(block_name)
        print("=" * len(block_name))
        for line in lines:
            for k, v in line.items():
                print(f"{k} = {v}")
        print()


# Script ============================================================================
if __name__ == "__main__":
    # Parse the script arguments
    parser = ArgumentParser(description="Parse arguments for main app")
    parser.add_argument(
        dest="conf_file",
        type=str,
        nargs="?",
        default=None,
        help="configuration file (namelist) for this script",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        type=str,
        required=False,
        default=".",
        help="select output directory for generated files",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose logging")
    parser.add_argument(
        "--template", action="store_true", help="print a configuration file template"
    )
    parser.add_argument(
        "--print_mcfost_defaults",
        action="store_true",
        help="print internal default values for all mcfost parameters as returned by vac2fost",
    )
    parser.add_argument("--version", action="store_true", help="display vac2fost version")

    cargs = parser.parse_args()

    if cargs.version:
        print(__version__)
        sys.exit(0)

    # special cases -----------------------------
    if cargs.template:
        print(template)
        print(f"! automatically generated with vac2fost {__version__}\n")
        sys.exit(0)
    elif cargs.print_mcfost_defaults:
        print_mcfost_default_conf()
        sys.exit(0)
    elif len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # -------------------------------------------
    main(
        cargs.conf_file,
        output_dir=cargs.output_dir.strip(),
        loglevel={True: 0, False: 30}[cargs.verbose],
    )
