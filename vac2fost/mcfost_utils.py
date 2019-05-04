"""Utility functions to call MCFOST in vac2fost.main()
to define the output grid."""
import os
import shutil
from pathlib import Path
from subprocess import run, CalledProcessError
from collections import OrderedDict as od
from uuid import uuid4 as uuid

import numpy as np
from astropy import fits


MINGRAINSIZE_µ = 0.1

blocks_descriptors = od(
    # This nested orderded dictionnary describes a default parafile for mcfost
    #
    # parameter names should match mcfost"s documentation
    # http://ipag-old.osug.fr/~pintec/mcfost/docs/html/parameter_file.html
    #
    # this is still WIP (remaining naming discrepencies)
    #
    # notes for a future pull-request on mcfost documentation itself:
    # *   indicates a change in name for various reasons...
    # ?   indicates a missing documentation line (or a deprecated parameter)
    # $   indicates stuff I"ll have to go over again, either because it
    #     breaks regression here, or because I need to change to api altogether

    # fixplan : DONE 1) PR to fix <> and % in different commits (and %% ?)
    #                2) ask Christophe about "?" and go over *
    #                3) deal with *
    #                4) deal with $

    [
        ("Photons", (
            od([("nbr_photons_eq_temp", "1.28e5")]),
            od([("nbr_photons_lambda", "1.28e3")]),
            od([("nbr_photons_image", "1.28e5")])
        )),
        ("Wavelengths", (
            od([("n_lambda", 50),
                ("lambda_min", "0.1"),
                ("lambda_max", "3e3")]),
            od([("ltemp", True),
                ("lsed", True),
                ("use_default_wavelength_grid", True)]),
            od([("wavelength_file", "wavelengths.dat")]),
            od([("separate_contributions", False),
                ("output_stokes_parameters", False)])
        )),
        ("Grid", (
            od([("grid_type", 1)]),
            od([("n_rad", 100),
                ("nz", 10),
                ("n_az", 100),
                ("n_rad_in", 30)])
        )),
        ("Images", (
            od([("grid_nx", 501),
                ("grid_ny", 501),
                ("map_size", 400)]),
            od([("RT_imin", 0),
                ("RT_imax", 0),
                ("RT_n_incl", 1),
                ("RT_centered", False)]),
            od([("RT_az_min", 0),
                ("RT_az_max", 240),
                ("RT_n_az", 1)]),
            od([("distance", 140)]),
            od([("disk_PA", 0)])
        )),
        ("Scattering Method", (
            od([("scattering_method", 0)]),
            od([("Mie_hg", 1)]) # *
        )),
        ("Symmetries", (
            od([("image_symmetry", True)]),
            od([("central_symmetry", True)]),
            od([("plane_symmetry", True)]),
        )),
        ("Disk physics", (
            od([("dust_settling", 0),
                ("exp_strat", 0.5),
                ("a_srat", 1.0)]),
            od([("dust_radial_migration", False)]),
            od([("sublimate_dust", False)]), # mcfost issue : check order !
            od([("hydrostatic_equilibrium", False)]),
            od([("viscous_heating", False),
                ("viscosity", "1e-3")]),
        )),
        ("Number of Zones", (
            od([("n_zones", 1)]), # *
        )),
        ("Density structure", (
            od([("zone_type", 1)]),
            od([("disk_dust_mass", "1e-3"),
                ("gas_to_dust_ratio", 100)]),
            od([("scale_height", 5.0), # *
                ("reference_radius", 100.0), # *
                ("vertical_profile_exponent", 2)]),
            od([("rin", 10), # $
                ("edge", 0),
                ("rout", 200), # $
                ("Rc", 100)]),
            od([("flaring_exp", 1.0)]), # *
            od([("density_exp", -0.5), # *
                ("gamma_exp", 0.0)]) # *
        )),
        ("Grain properties", (
            od([("n_species", 1)]),
            od([("Grain_type", "Mie"),
                ("n_components", 1),
                ("mixing_rule", 2),
                ("porosity", 0.),
                ("mass_fraction", 0.75),
                ("DHS_Vmax", 0.9)]),
            od([("optical_indices_file", "Draine_Si_sUV.dat"),
                ("volume_fraction", 1.0)]),
            od([("heating_method", 1)]),
            od([("amin", MINGRAINSIZE_µ),
                ("amax", 1000),
                ("aexp", 3.5),
                ("n_grains", 100)])
        )),
        ("Molecular RT settings", (
            od([("lpop", True),
                ("lpop_accurate", True),
                ("LTE", True),
                ("profile_width", 15.0)]),
            od([("v_turb", 0.0)]), # ?
            od([("nmol", 1)]), # ?
            od([("molecular_data_file", "13co.dat"), # *
                ("level_max", 6)]),
            od([("vmax", 1.0),
                ("n_speed", 20)]),
            od([("cst_abundance", True),
                ("abund", "1e-6"), # ?
                ("abund_file", "abundance.fits.gz")]), # ?
            od([("ray_tracing", True), # *
                ("n_lines", 3)]),
            od([("transition_num_1", 1),
                ("transition_num_2", 2),
                ("transition_num_3", 3)])
        )),
        ("Star properties", (
            od([("n_stars", 1)]),
            od([("Teff", 4000.0),
                ("Rstar", 2.0),
                ("Mstar", 1.0),
                ("x", 0.),
                ("y", 0.),
                ("z", 0),
                ("is_blackbody", True)]),
            od([("star_rad_file", "lte4000-3.5.NextGen.fits.gz")]), # ?
            od([("fUV", 0.0), ("slope_fUV", 2.2)]),
        ))
    ])

KNOWN_MCFOST_ARGS = []
for descriptor in blocks_descriptors.items():
    for di in descriptor[1]:
        KNOWN_MCFOST_ARGS += [k.lower() for k in di.keys()]

def write_mcfost_conf(output_file: Path, custom: dict = None, verbose=False):
    """Write a configuration file for mcfost using values from <custom>,
    and falling back to defaults found in block_descriptor defined above
    """
    if custom is None:
        custom = {}
    if output_file.exists() and verbose:
        print(f"Warning: {output_file} already exists, and will be overwritten.")
    with open(output_file, mode="wt") as fi:
        fi.write(".".join(min_mcfost_version.split(".")[:2]).ljust(10) +
                 "mcfost minimal version prescribed by vac2fost\n\n")
        for block, lines in blocks_descriptors.items():
            fi.write(f"# {block}\n")
            for line in lines:
                parameters = []
                for param, default in line.items():
                    if param in custom:
                        val = custom[param]
                    else:
                        val = default
                    parameters.append(str(val))
                fi.write("  " + "  ".join(parameters).ljust(36)
                         + "  " + ", ".join(line.keys()))
                fi.write("\n")
            fi.write("\n")
        fi.write("\n\n")
        fi.write(f"%% automatically generated with vac2fost {__version__}\n")
        fi.write(f"%% via mcfost {DETECTED_MCFOST_VERSION}\n")
        fi.write(f"%% run by {os.environ['USER']} on {gethostname()}\n")
    if verbose:
        print(f"wrote {output_file}")

def get_mcfost_grid(itf) -> np.ndarray:
    """Pre-run MCFOST with -disk_struct flag to get the exact grid used."""
    mcfost_conf_file = itf.mcfost_conf_file
    output_dir = itf.io.OUT.directory

    output_dir = Path(output_dir).resolve()
    mcfost_conf_path = Path(mcfost_conf_file)
    if not output_dir.exists():
        os.makedirs(output_dir)

    grid_file_name = output_dir / "mcfost_grid.fits.gz"

    if itf.current_num == itf.nums[0]:
        assert mcfost_conf_path.exists()
        # generate a grid data file with mcfost itself and extract it
        tmp_mcfost_dir = output_dir / f"TMP_VAC2FOST_MCFOST_GRID_{uuid()}"
        os.makedirs(tmp_mcfost_dir)
        try:
            shutil.copyfile(mcfost_conf_path.resolve(),
                            tmp_mcfost_dir/mcfost_conf_path.name)
        except shutil.SameFileError:
            pass

        pile = Path.cwd()
        os.chdir(tmp_mcfost_dir)
        try:
            run(["mcfost", "mcfost_conf.para", "-disk_struct"], check=True,
                capture_output=(not itf.mcfost_verbose))

            shutil.move("data_disk/grid.fits.gz", grid_file_name)
        except CalledProcessError as exc:
            errtip = f"\nError in MCFOST, exited with exitcode {exc.returncode}"
            if exc.returncode == 174:
                errtip += (
                    "\nThis is probably a memory issue. "
                    "Try reducing the target resolution or,"
                    " alternatively, give more cpu memory to this task."
                )
            raise RuntimeError(errtip) from exc
        finally:
            os.chdir(pile)
            shutil.rmtree(tmp_mcfost_dir)
    with fits.open(grid_file_name, mode="readonly") as fi:
        target_grid = fi[0].data
    return target_grid
