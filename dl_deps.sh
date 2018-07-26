#/usr/bin/sh

set -eu
if [ -d deps ] ; then
    rm -fr deps
fi

mkdir deps
cd deps

git clone https://gitlab.oca.eu/crobert/amrvac-pywrap-project.git
git clone https://gitlab.oca.eu/crobert/vtk_vacreader-project.git
