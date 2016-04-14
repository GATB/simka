#!/bin/bash
#simple test with real data

# look for simka binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/simka" ]
then
 bindir="../bin"
elif [ -f "../build/bin/simka" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled simka binary"
 exit 1
fi

# run simka
$bindir/simka -verbose 0 -in ../example/simka_input.txt -out output/ -out-tmp temp_output
var=$?
if [ $var -eq 0 ]
then
    echo  "*** Test: PASSED"
else
    echo  "*** Test: FAILED"
    exit 1
fi

# clean temp files
rm -rf output temp_output
