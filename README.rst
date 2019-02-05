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

The recommended method relies on the package management system and environment management system ``conda``.

!!!In some version of conda, the ``conda develop`` command ignore the environment python version. The root environment version must be >3. If you use modules on your computer/cluster, you should load ``python3``.

The following assumes ``currrent working directory == vac2fost-project``.

Most of the non-standard Python dependencies can be installed with

    .. code-block:: bash
    
        conda create --name vac2fost --file environment.yml --channel conda-forge

This will create a specific environment for using vac2fost without modifying your usual python packages in your standard environment.
The two sub-parts of the program are amrvac_pywrap_ and vtk_vacreader_.

.. _amrvac_pywrap: https://gitlab.oca.eu/crobert/amrvac-pywrap-project
.. _vtk_vacreader: https://gitlab.oca.eu/crobert/vtk_vacreader-project

They can be installed with

    .. code-block:: bash

        source activate vac2fost
        bash install_deps.sh # this script will ask for user confirmation

The first line change your working environment to the one needed to use vac2fost. The second will install the two sub-programs.

Then, there are two possibilities to use ``vac2fost``:

- if you wish to use ``vac2fost`` as a Python package, install it as

    .. code-block:: bash

        source activate vac2fost
        conda install .
        #or
        conda develop . # if you wish to actively modify the package as you use it


- if you wish to use ``vac2fost.py`` from command line, the recommended fashion is to create a symbolic link to the main file, as part of your ``$PATH`` and treat it as an executable. For instance: 

    .. code-block:: bash
        
        chmod +x vac2fost/vac2fost.py
        ln -s vac2fost/vac2fost.py ~/local/bin/vac2fost.py
        #or a more standard solution
        export PATH=$PATH:~/path/to/vac2fost



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

           &amrvac_input
            config = 'relative/to/<hydro_data_dir>/path/to/amrvac/config/file1.par','and/file2.par'
            conv2au = 100
            hydro_data_dir = 'path/to/output/data/directory'
            nums = 0
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


The app can be used in two fashions

* directly from command-line:

  .. code:: bash

            # provided that the num parameter is included in the configuration:&amrvac_input:nums
            ./vac2mcfost.py <configuration_file> --dbm <[dust-only, gas-only, mixed]>
            # otherwise
            ./vac2mcfost.py <configuration_file> --nums <input file num>

* as an importable python function

  .. code:: python

            from vac2fost import main as vac2fost

            conf = ...  #(str or pathlib.Path)
            out = ...   #(str or pathlib.Path)

	    # minimal call
            vac2fost(config_file=conf)

	    # more sophisticated call
            vac2fost(config_file=conf, nums=10, output_dir=out)
  
note that if ``nums`` are defined as a parameter **and** included in
the configuration, the parameter value is used.
``nums`` can be a single integer or any integer-returning iterable.

Get help
--------

To see optional parameters available, run

  .. code:: bash

	    vac2fost.py --help
