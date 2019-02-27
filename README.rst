.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat

warning: this work has seen very little usage outside of its author's
laptop, any user feedback is therefor essential to making this program
more reliable. Please do report bugs on sight, thanks.

vac2fost
========

vac2fost is a Python 3.6+ program that translates `.vtu`
formated MPI-AMRVAC output files into `.fits` files compatible with
mcfost 3D model input format.

``vac2fost.py`` can be used for python scripting, importing its ``main()``
function within Python, or as an command-line executable.


Content
-------

``tests/sample/`` contains sample configuration files used by test
scripts located in ``tests/``

Any mcfost parameter can be overwritten by those found in
``&mcfost_list``, using names defined in
``vac2fost.py:MCFOSTUtils:block_descriptors``.


Prerequisites
-------------

As vac2fost calls mcfost itself in order to define the output grid, it
is assumed that you have ``mcfost (>=3.0.35)`` properly install.

The recommended method for installing vac2fost consists in creating a
dedicated and _isolated_ Python environment.  `environment.yml`
describes dependencies and the rationale for each version constraint.

This file can be used to create a environment in one line with conda

    .. code-block:: bash
    
        conda create --name vac2fost --file environment.yml --channel conda-forge

For details, see the official conda documentation for managing environments
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

The program also relies on a derived project, vtk_vacreader_: a small
module to read vtk files from MPI-AMRVAC (restricted to non-AMR files
atm).

.. _vtk_vacreader: https://gitlab.oca.eu/crobert/vtk_vacreader-project

It can can be installed as

    .. code-block:: bash

        conda activate vac2fost
	# old versions of anaconda use "source" instead of "conda" in this command

        bash install_deps.sh # this script will ask for user confirmation

**warning; old versions of conda contain a bug where `conda develop` breaks if running in a environement with a different Python version than the base one. This is solved by updating conda with `conda update conda``
If for whatever reason you are unable to do that, vtk_vacreader can be installed with pip, as long as you're using an environement specific pip.
**



Usage/Finishing installation
----------------------------

vac2fost can be used in two fashions:

- as a Python package. Install with

    .. code-block:: bash

        conda activate vac2fost
        conda install .
        #or
        conda develop . # if you wish to actively modify the package as you use it


- from command line. The recommended fashion is to create a symbolic
  link to the main file, as part of your ``$PATH`` and treat it as an
  executable.  For instance

    .. code-block:: bash
        
        chmod +x vac2fost/vac2fost.py
        ln -s vac2fost/vac2fost.py ~/local/bin/vac2fost.py
        #or a more standard solution
        export PATH=$PATH:~/path/to/vac2fost



Testing
-------

run ``pytest`` to check that the program behaves normally. Tests
should be able to pass as long as your environment contains all of vac2fost's dependencies,
even if you did not install vac2fost`` yet.

A demo of how to produce an image with ``vac2fost`` + ``MCFOST`` can
be found in ``docs/mcfost_demo.sh``.


Getting started
---------------

The minimal requirement is a configuration file, which 
can be generated from command line with ``vac2fost.py --genconf > conf.nml``

For instance
 .. code:: fortran

	   &amrvac_input
	   config = 'relative/to/<hydro_data_dir>/path/to/amrvac/config/file1.par','and/file2.par'
           hydro_data_dir = 'path/to/output/data/directory'
           nums = 0
           /

	   &units
	   ! conversion factors between dimensionel amrvac outputs and physical units
	   distance2au = 100.0
	   time2yr     = 10.
	   /

	   &mcfost_output
           ! this list describes MCFOST parameters
           ! named according to vac2fost.MCFOSTUtils.blocks_descriptors
           nr   = 150
           nphi = 100
           nz   = 50
           nr_in = 30  ! need to be < nr

           flaring_index = 1.125
           ref_radius = 100.0    ! [a.u.]
           scale_height = 10.0   ! [a.u.] defined at ref_radius

           star_mass = 1.8
           star_temp = 6550
           distance  = 157
	   /

How to use it

* from command-line:

  A typical call would look like this
  .. code:: bash

            # provided that the num parameter is included in the configuration:&amrvac_input:nums
            ./vac2mcfost.py <configuration_file> --dbm <[dust-only, gas-only, mixed]>
            # otherwise
            ./vac2mcfost.py <configuration_file> --nums <input file num(s)>

* as an importable python function

  .. code:: python

            from vac2fost import main as vac2fost

            conf = ...  #(str or pathlib.Path)
            out = ...   #(str or pathlib.Path)

	    # minimal call
            vac2fost(config_file=conf)

	    # more sophisticated call
            vac2fost(config_file=conf, nums=10, output_dir=out)
  
note that if ``nums`` are defined as a command line arguemnt **and**
included in the configuration file, the argument prevails.  ``nums``
can be a single integer or any integer-returning iterable.

Dust binning mode
-----------------

va2fost can be used wether or not your hydro simulation contains dust.
The way it works is by guessing the most appropriate thing to do,
encoded in a parameter called `dust-binning-mode` (or "dbm" for
shorts):

- if no dusty fluid is found, gas will be used as a proxy, and all
  grain sizes will be assumed to follow gas distribution
  (`dbm="gas-only")
- if dust is found but no one species is smaller than 0.1 micron, gas
  is still used to trace the smallest grains (`dbm="mixed")

By default, vac2fost automatically sets the dbm, but it can be imposed
by the user as an argument.  An additional mode is "dust-only", where
gas density is being ignored. This mode is never chosen automatically
but can prove relevant for tests.

If dbm is set to "dust-only", one can also pass gas density as gas
itself to mcfost with "read_gas_density". In other dbms, this
parameter is ignored because mcfost is already assuming that gas and
smallest grains are perfectly coupled.


Get help
--------

vac2fost's command line help is displayed upon
  .. code:: bash

	    vac2fost.py --help

	    #or even simplier
	    vac2fost.py

