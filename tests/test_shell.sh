#! /usr/bin/sh

python -c "import vac2fost"
if [ $? == 1 ] ; then
    echo "error: vac2fost is not installed properly"
    exit 1
fi

if [ -z ${VAC2FOST_INSTALL_DIR+x} ] ; then
    echo "error: this script requires VAC2FOST_INSTALL_DIR to be defined."
    exit 1
fi
echo "install dir used : $VAC2FOST_INSTALL_DIR"

TESTDIR=$VAC2FOST_INSTALL_DIR/tests
OUT=$TESTDIR/output
EXE=$VAC2FOST_INSTALL_DIR/app/v2f.py

expect_success () {
    if [[ $? != 0 ]] ; then
	echo FAIL
        exit 1
    fi
    echo "PASS"
}

expect_failure () {
    if [[ $? == 0 ]] ; then
	echo FAIL
        exit 1
    fi
    echo "PASS (expected failure)"
}


mkdir -p $OUT

echo "Test 1: normal call"
$EXE $TESTDIR/sample/vac2fost_conf_quick.nml --output $OUT/shell_1
expect_success

# <- this one is SUPPOSED to fail because the mandatory "nums" argument is not provided
echo
echo "Test 2: Broken call"
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_2
expect_failure

echo "Test 3: display version"
$EXE --version
expect_success
