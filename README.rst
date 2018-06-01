VAC2FOST
========

`vac2mcfost.py` is a `python3` interface script that translates a `.vtu`
amrvac output file into a `.fits` that can be fed to mcfost through
the option `-density_file`.


.. note::

   What I still need to do here:

   * add more tests
        a) handle lists of conf (already supported but currently not being tested)
        b) call mcfost
   * consider adding switch for cylindrical/spherial griding
        a) handle at least one case correctly (with zmax and scale_height issues)
   * add option for dust settling method (Gauss vs Fromang, cf McFost options)



Dependencies
------------

This package relies on `amrvac_pywrap` and `vtk_vacreader` packages (developed
by Clément Robert).



Content
-------



`tests/sample/` contains sample configuration files used by `tests/test_app.py` 
and `tests/test_shell.py` 

`vac2fost/data/default_mcfost_conf.para` contains a basic mcfost configuration
from which is generated `mcfost_conf.para` (into the current directory) when the
app is called, where certain parameters can be overwritten with value found in
`api_script.nml:&mcfost_list`



Usage
-----

Run `pytest` in this directory to check if the interface works in your 
environment, then `./gen_image.sh` to generate an actual image with mcfost.

The minimal requirement is a `configuration.nml` (name does not have
to match this example) file formated as follow

 .. code:: fortran

           &mcfost_list
           ! this list describes mcfost parameters
           ! only a clearly identified set of them can be passed this way
           ! those are found in the python dictionnary mcfost_args_locate in vac2fost.py
               nr   = 150
               nphi = 100
               nz   = 50
               nr_in = 30  ! need to be < nr
               star_mass = 1.8  ! [M_sun] (Casassus et al 2015)
               star_temp = 6550 ! [K]     (Mendiguita et al 2014)
               distance  = 157  ! [pc]    (Gaia consortium 2016)
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
            ./vac2mcfost.py $configuration_file
            # otherwise
            ./vac2mcfost.py $configuration_file -o $offset

* as an importable python function

  .. code:: python

            from vac2fost import main as vac2fost

            conf_file = ... #(str or pathlib.Path)
            dat_file  = ... #(str or pathlib.Path)
            
            vac2fost(conf_file)
            vac2fost(conf_file, offset)
  
note that if `offset` is defined as a parameter **and** included in
the configuration, the parameter value is used.
