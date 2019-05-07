"""Shell level script"""
import sys
from pathlib import Path
from argparse import ArgumentParser

import f90nml

from vac2fost.info import __version__
from vac2fost.utils import CYAN, BOLD, get_prompt_size, decorated_centered_message
from vac2fost.interfaces import DEFAULT_UNITS, Interface, VerbatimInterface




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


def main(config_file: str,
         nums: int = None, # or any in-returning interable
         output_dir: Path = Path.cwd(),
         dust_bin_mode: str = "auto",
         read_gas_density=False,
         read_gas_velocity=False,
         verbose=False,
         mcfost_verbose=False):
    """Transform a .vtu datfile into a .fits"""

    print(decorated_centered_message("start vac2fost"))

    InterfaceType = {True: VerbatimInterface, False: Interface}[verbose]
    itf = InterfaceType(config_file, nums=nums, output_dir=output_dir,
                        dust_bin_mode=dust_bin_mode,
                        read_gas_density=read_gas_density,
                        read_gas_velocity=read_gas_velocity,
                        mcfost_verbose=mcfost_verbose)

    for i, n in enumerate(itf.nums):
        if verbose or i == 0:
            print()
        mess1 = f"current input number: {n}"
        mess2 = f"({i+1}/{len(itf.nums)})"
        print((BOLD+" "*(get_prompt_size()-len(mess1)-len(mess2))).join([mess1, mess2]))
        print("-"*get_prompt_size())
        try:
            itf.load_input_data(n)
        except FileNotFoundError as err:
            filepath = Path(str(err)).relative_to(Path.cwd())
            itf.warnings.append(f"file not found: {filepath}")
            continue
        itf.write_mcfost_conf_file()
        itf.gen_2D_arrays()
        itf.gen_3D_arrays()
        itf.write_output()

        try:
            filepath = itf.io.OUT.filepath.relative_to(Path.cwd())
        except ValueError:
            filepath = itf.io.OUT.filepath
        print(CYAN + f" >>> wrote {filepath}")

    if itf.warnings:
        print()
        itf.display_warnings()

    print(decorated_centered_message("end vac2fost"))

    # return the Interface object for inspection (tests)
    return itf

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
        "--read_gas_velocity",
        action="store_true",
        help="pass gas velocity to mcfost (keplerian velocity is assumed otherwise)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="activate verbose mode"
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

    if cargs.genconf:
        print(generate_conf_template())
        print(f"%% automatically generated with vac2fost {__version__}\n")
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
        verbose=cargs.verbose,
        mcfost_verbose=cargs.mcfost_verbose
    )
    # -------------------------------------------
    if cargs.cprofile:
        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
        ps.print_stats()
        print(s.getvalue())
