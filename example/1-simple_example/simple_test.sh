#!/bin/bash
#simple test with real data

# look for simka binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
#if [ -f "../../scripts/simka2/bin/simka" ]
#then
# bindir="../bin"
#elif [ -f "../build/bin/simka" ]
#then
# bindir="../build/bin"
#else
# echo "could not find a compiled simka binary"
# exit 1
#fi

# run simka
command="python ../scripts/simka2/simka.py -in ../example/1-simple_example/simka_input.txt -out output/ -out-tmp temp_output"
#printf "$command\n\n"

$command

printf "\n\n\n"
var=$?
if [ $var -eq 0 ]
then
    echo  "*** Test: PASSED"
else
    echo  "*** Test: FAILED"
    exit 1
fi

#printf "\nremoving all created dirs\n"
# clean temp files
rm -rf temp_output

printf "\nCommand used:\n" 
printf "\t$command\n"
printf "Command for creating heatmaps and dendrogram:\n"
printf "\tpython ../scripts/create_heatmaps.py output/\n"
