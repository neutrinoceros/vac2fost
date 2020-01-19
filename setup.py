from setuptools import setup, find_packages
from pathlib import Path
import sys

if sys.version < "3.7":
    raise OSError("vac2fost is written in Python 3.7")
elif sys.version < "3.7.1":
    raise OSError(
        "Python 3.7.0 has known issues with scipy, please update to Python 3.7.1 or later"
    )

here = Path(__file__).parent.absolute()
project_name = "vac2fost"

setup(
    name=project_name,
    version=__import__(project_name).__version__,
    author=__import__(project_name).__author__,
    author_email=__import__(project_name).__contact__,
    url="https://gitlab.oca.eu/crobert/vac2fost-project",
    classifiers=["Development Status :: 3 - Alpha", "Programming Language :: Python"],
    keywords="interface data-analysis",
    install_requires=[
        "scipy>=1.1",
        "numpy>=1.15",
        "astropy>=3.1.1",
        "vtk>=8.1.2",
        "f90nml>=1.0.2",
    ],
    python_requires=">=3.7",
    packages=find_packages(),
    package_data={"vac2fost": ["data/default_mcfost_conf.para"]},
)
