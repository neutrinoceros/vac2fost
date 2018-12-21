#!/usr/bin/env bash

# vac2fost demo script.
# ======================================================
# generate an image with mcfost from AMRVAC
# density data .vtu input in 3 steps
#   1. preroll : run vac2fost.py
#   2. run MCFOST to compute temperature
#   3. run MCFOST to calculate an image at 1300 microns
# ======================================================

# globals ------------------------------------------------------------------
OUT=demo_out
DENS=$OUT/hd142527_dusty0000.fits
CONF=$OUT/mcfost_conf.para
MCFOST_OUT=$OUT/mcfost_out

# preliminary checks -------------------------------------------------------
MY_MCFOST=$(which mcfost)
if [ $? == 0 ] ; then
    echo "Found mcfost installed under $MY_MCFOST"
else
    echo "Could not locate mcfost installation."
    exit 1
fi

MY_APP=$(which vac2fost.py)
if [ $? == 0 ] ; then
    echo "Found vac2fost installed under $MY_APP"
elif [ -f vac2fost/vac2fost.py ] ; then
    MY_APP=vac2fost/vac2fost.py
    echo "Found vac2fost installed under $MY_APP"
else
    echo "Could not locate mcfost installation."
    exit 2
fi

# preroll the interface ----------------------------------------------------
$MY_APP tests/sample/vac2fost_conf.nml -o $OUT --verbose

# inter-code checks --------------------------------------------------------
[ ! -d "$MCFOST_OUT" ] && mkdir $MCFOST_OUT

if [ ! -f "$DENS" ]; then
    echo "input file $DENS not found, aborting demo"
    exit 3
fi

if [ ! -f "$CONF" ]; then
    echo "$CONF not found, aborting demo"
    exit 4
fi

# compute temperature ------------------------------------------------------
MCFOST_CALL="mcfost $CONF -density_file $DENS -root_dir $MCFOST_OUT -3D"
eval $MCFOST_CALL

if [ $? != 0 ] ; then
    echo "first call to mcfost failed, aborting demo"
    exit 5
fi

# compute image ------------------------------------------------------------
eval $MCFOST_CALL -img 1300
