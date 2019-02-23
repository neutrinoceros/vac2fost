#! /usr/bin/sh
TESTDIR=$ROOTDIR/tests
OUT=$TESTDIR/output
EXE=$OUT/v2f_exe.py

mkdir -p $OUT
echo $EXE
cp $ROOTDIR/vac2fost/vac2fost.py $EXE
chmod +x $EXE


$EXE $TESTDIR/sample/vac2fost_conf_quick.nml --output $OUT/shell_1
$EXE $TESTDIR/sample/vac2fost_conf_quick.nml --output $OUT/shell_2 --nums 2
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_3 --nums 0
$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_4 --nums 0 1 2
#$EXE $TESTDIR/sample/vac2fost_conf_quick_no_number.nml --output $OUT/shell_5 # <- this one is SUPPOSED to fail

