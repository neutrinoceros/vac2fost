VAC2FOST
========

``vac2mcfost.py`` is a ``Python 3`` interface script that translates a
``.vtu`` ``MPI-AMRVAC`` output file into a ``.fits`` that can be fed to
``MCFOST`` through the option ``-density_file``.


Python uncommon dependencies
----------------------------

This package relies on ``amrvac_pywrap`` and ``vtk_vacreader`` Python 3
packages (developed by Clément Robert).


Other dependencies
------------------

Python:

    - astropy 3.1.1 (to come)
    - vtk 8 (only compatible with python 3.6, will be dropped at some point)
    - f90nml

Misc:

    - ``MCFOST 3`` is used soly to compute the target grid. Hopefully we
      can drop that dependency in the future


Content
-------

``tests/sample/`` contains sample configuration files used by test
scripts located in ``tests/``

``vac2fost/data/default_mcfost_conf.para`` contains a basic mcfost
configuration from which ``mcfost_conf.para`` is generated in the output
directory when ``vac2fost.main()`` is called.  Any ``MCFOST`` parameter
can be overwritten by those found in ``&mcfost_list``, using names
defined in ``vac2fost.py``.


Installation
------------

The program can be used either as a Python package with a ``main()``
fonction, or from command line.

To install it as a package within anaconda, use ``conda develop
<path/to/vac2fost-project/>``.

To use the program as a script from command line, the main file should
be accessible through $PATH.  A relatively clean way to do this is by
linking a dummy file to the installation directory.


Testing
-------

run ``pytest`` to check that the script behaves normally. Tests should
pass even if you did not install ``vac2fost`` yet, either as a package
or an executable.

A demo of how to produce an image with ``vac2fost`` + ``MCFOST`` can
be found in ``docs/mcfost_demo.sh``.


Usage
-----

The minimal requirement is a ``configuration.nml`` (name does not have
to match this example) file formated as follow (such a file can be
generated from command line with ``vac2fost.py --genconf``)

 .. code:: fortran

           &mcfost_list
           ! this list describes MCFOST parameters
           ! named according to vac2fost.MCFOSTUtils.blocks_descriptors
               nr   = 150
               nphi = 100
               nz   = 50
               nr_in = 30  ! need to be < nr
               star_mass = 1.8
               star_temp = 6550
               distance  = 157
           /

           &target_options
           ! additional options
               origin = '/path/to/mod_usr.t/parent/directory'
               amrvac_conf = 'relative/path/to/vac/config_file/from/origin'
               offset = 0  ! output number of the .dat file to be converted
               zmax = 5    ! use same unit at distance in the original simulation
               aspect_ratio = 0.01
           /


The app can be used in two fashions

* directly from command-line:

  .. code:: bash

            # provided that the offset is included in the configuration:&target_options:offset
            ./vac2mcfost.py <configuration_file>
            # otherwise
            ./vac2mcfost.py <configuration_file> -n <input file num>

* as an importable python function

  .. code:: python

            from vac2fost import main as vac2fost

            conf = ...  #(str or pathlib.Path)
            out = ...   #(str or pathlib.Path)

	    # minimal call
            vac2fost(config_file=conf)

	    # more sophisticated call
            vac2fost(config_file=conf, offset=10, output_dir=out)
  
note that if ``<input file num>`` (aka ``offset``) is defined as a
parameter **and** included in the configuration, the parameter value
is used.


Get help
--------

To see optional parameters available, run

  .. code:: bash

	    vac2fost.py --help
