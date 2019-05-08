#! /usr/bin/sh

python -c "import vac2fost"
if [ $? == 1 ] ; then
    echo "error: vac2fost ins not installed properly"
    exit 1
fi

if [ -z ${VAC2FOST_INSTALL_DIR+x} ] ; then
    echo "error: this script requires VAC2FOST_INSTALL_DIR to be defined."
    exit 1
fi
echo "install dir used : $VAC2FOST_INSTALL_DIR"

TESTDIR=$VAC2FOST_INSTALL_DIR/tests
OUT=$TESTDIR/output
EXE=$OUT/v2f_exe.py


expect_success () {
    if [[ $? != 0 ]] ; then
	echo FAIL
        exit 1
    fi
    echo "PASS"
}

expect_faillure () {
    if [[ $? == 0 ]] ; then
	echo FAIL
        exit 1
    fi
    echo "PASS (expected faillure)"
}


mkdir -p $OUT
echo $EXE
cp $VAC2FOST_INSTALL_DIR/vac2fost/vac2fost.py $EXE
chmod +x $EXE

echo "Test 1: normal call"
$EXE $TESTDIR/sample/vac2fost_conf_quick.nml --output $OUT/shell_1
expect_success

echo
echo "Test 2: nums arg"
$EXE $TESTDIR/sample/vac2fost_conf_quick.nml --output $OUT/shell_2 --nums 2
expect_success

echo
echo "Test 3: nums = 0"
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_3 --nums 0
expect_success

echo
echo "Test 4: multinums call"
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_4 --nums 0 1 2
expect_success

# <- this one is SUPPOSED to fail
echo
echo "Test 5: Broken call"
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_5
expect_faillure
