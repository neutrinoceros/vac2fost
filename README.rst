VAC2FOST
========

``vac2fost.py`` is a ``Python 3`` interface script that translates a
``.vtu`` ``MPI-AMRVAC`` output file into a ``.fits`` that can be fed to
``MCFOST`` through the option ``-density_file``.


``vac2fost.py`` can be used either as a Python package with a ``main()``
fonction, or from command line.


Content
-------

``tests/sample/`` contains sample configuration files used by test
scripts located in ``tests/``

``vac2fost/data/default_mcfost_conf.para`` contains a basic mcfost
configuration from which ``mcfost_conf.para`` is generated in the output
directory when ``vac2fost.main()`` is called.

Any ``MCFOST`` parameter can be overwritten by those found in ``&mcfost_list``,
using names defined in ``vac2fost.py``.


Installation
------------
It is currently required to build a very specific Python environment in order to
run the program.

It is also assumed that you have ``mcfost (>=3.0)`` properly install

The recommended method relies on ``conda``.

The following assumes ``cwd == vac2fost-project``.

Most of the non-standard Python dependencies can be installed with

    .. code-block:: bash
    
        conda config --add channels conda-forge
        conda create --name vac2fost --file environment.ylm

Other parts of the program are ``amrvac_pywrap`` and ``vtk_vacreader``.

They can be installed with

    .. code-block:: bash

        conda activate vac2fost
        bash install_deps.sh # this script will ask for user confirmation


Then, if you wish to use ``vac2fost`` as a Python package, install it as

    .. code-block:: bash

        conda activate vac2fost
        conda install .
        #or
        conda develop . # if you wish to actively modify the package as you use it


Finally, if you wish to use ``vac2fost.py`` from command line, the recommended
fashion is to create a symbolic link to the main file, as part of your ``$PATH``
and treat it as an executable.

For instance: 

    .. code-block:: bash
        
        chmod +x vac2fost/vac2fost.py
        ln -s vac2fost/vac2fost.py ~/local/bin/vac2fost.py




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
               num = 0  ! output number of the .dat file to be converted
               zmax = 5    ! use same unit at distance in the original simulation
               aspect_ratio = 0.01
           /


The app can be used in two fashions

* directly from command-line:

  .. code:: bash

            # provided that the num parameter is included in the configuration:&target_options:num
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
            vac2fost(config_file=conf, num=10, output_dir=out)
  
note that if ``<input file num>`` is defined as a parameter **and** included in
the configuration, the parameter value is used.


Get help
--------

To see optional parameters available, run

  .. code:: bash

	    vac2fost.py --help
