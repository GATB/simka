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
command="$bindir/simka -in ../example/simka_input.txt -out ./simka_results/ -out-tmp ./simka_temp_output"
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

printf "\nCommand for visualizing results:\n"
printf "\tpython ../scripts/visualization/run-visualization.py -in ./simka_results/ -out ./simka_results/ -pca -heatmap -tree\n"

printf "\nCommand for visualizing results with metadata annotations:\n"
printf "\tpython ../scripts/visualization/run-visualization.py -in ./simka_results/ -out ./simka_results/ -pca -heatmap -tree -metadata-in ../example/dataset_metadata.csv -metadata-variable VARIABLE_1\n"