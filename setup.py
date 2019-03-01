from setuptools import setup, find_packages
from pathlib import Path
import os

here = Path(__file__).parent.absolute()
project_name = "vac2fost"

setup(
    name = project_name,
    version        = __import__(project_name).__version__,
    author         = __import__(project_name).__author__,
    author_email   = __import__(project_name).__contact__,
    url = "https://gitlab.oca.eu/crobert/vac2fost-project",
    license = "GNU",
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python"
    ],
    keywords = "interface data-analysis",
    install_requires = [
        "scipy>=1.1",
        "numpy>=1.15",
        "astropy>=3.1.1"
        "vtk>=8.1.2",
        #"f90nml>=1.0.2",
        #"vtk_vacreader>=1.0.1"
    ],
    python_requires = ">=3.7",
    packages = find_packages(),
    package_data = {"vac2fost": ["data/default_mcfost_conf.para"]},
)
