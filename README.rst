.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat

python macro
============

``vac2fost`` is a Python 3.7 macro that converts `.vtu` formatted
MPI-AMRVAC output files into `.fits` files compatible with mcfost 3D
model input format.

Installation (via conda)
------------------------

It is recommended to create a separate conda enviromnent before installing.
See https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

**prerequisite:** ``vac2fost`` currently relies on a derived project, vtk_vacreader_: a small
module to read VTK files from MPI-AMRVAC (restricted to non-AMR files
atm).
It can can be installed with the interactive script `install_deps.sh`.
Please make sure to run it from your dedicated python environment, be
it conda or pipenv... This script will ask for your prefered installer (pip or conda ?).


The minimal working installation script is

.. code-block:: bash

    bash install_deps.sh
    conda install --file environment.yml --channel conda-forge
    conda install .

It is recommended to replace the last line with ``conda develop . `` in case
you wish to actively modify the code while using it.



Usage
-----

.. code :: python

    from vac2fost import main as vac2fost
    # minimal requirements (will output to cwd)
    vac2fost("conf.nml")
    
    # specifying the output directory
    vac2fost("conf.nml", output_dir="my/out/")
    
    # overriding or adding some properties to "conf.nml"
    vac2fost("conf.nml", override={"flags": {"axisymmetry": True}})



Testing
-------

Run ``pytest`` from the installation directory to check that the program behaves
normally.




Complementary command line tool
===============================

Installation
------------
 **note** the command line tool relies on the package being installed
 
The file is found in the ``app/`` folder.
The recommended installation method is to alias it with a symbolic link such as:

    .. code-block:: bash
        ln -s path/to/vac2fost-project/app/v2f.py /usr/bin/vac2fost

Usage
-----

.. code:: shell

    # echo a conf file template
    vac2fost --template
    
    # display mcfost default parameters (names and values)
    vac2fost --print_mcfost_defaults
    
    # show version, help
    vac2fost --version
    vac2fost --help





How vac2fost selects density fields
===================================

``va2fost`` can be used wether or not your hydro simulation contains dust.
The way it works is by guessing the most appropriate thing to do,
encoded in a parameter called `dust-binning-mode` (or *DBM* for
shorts):

- if ``&dust`` namelist is ommited in the configuration file,
  gas will be used as a proxy, and all grain sizes will be assumed to follow gas
  distribution (DBM = "gas-only")
- if dust is found but no one species is smaller than 0.1 micron, gas
  is still used to trace the smallest grains (DBM = "mixed")

By default, vac2fost automatically sets the dbm, but it can be overrided via
``&flags: dust_bin_mode``. An additional mode is "dust-only", where
gas density is being ignored. This mode is never chosen automatically
but can prove relevant for tests.

If dbm is set to "dust-only", one can also pass gas density as gas
itself to mcfost with "read_gas_density". Within other DBMs, this
parameter is ignored because mcfost is already assuming that gas and
smallest grains are perfectly coupled.


.. _vtk_vacreader: https://gitlab.oca.eu/crobert/vtk_vacreader-project