#/usr/bin/sh

set -eu
if [ -d deps ] ; then
    read -p "Found existing deps directory. Overwrite it ? [y]/n    " choice
    case $choice in
	[nN]*)
	    echo "exiting install script !"
	    exit 1
	    ;;
	*)
	    rm -fr deps
    esac
fi

mkdir deps
cd deps
git clone https://gitlab.oca.eu/crobert/amrvac-pywrap-project.git
git clone https://gitlab.oca.eu/crobert/vtk_vacreader-project.git

read -p "Do you wish to install deps within current conda env ? y/[n]    " choice

case $choice in
    [yY]*) 
	conda develop amrvac-pywrap-project
	conda develop vtk_vacreader-project
	;;
    *)
	echo "deps were downloaded but not installed !"
esac

