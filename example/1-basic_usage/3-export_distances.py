

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])





export_bin = "../../scripts/simka2/bin/simka2-export"
input_dir = "simka_results/matrix_binary"

if not os.path.exists(input_dir):
    print("Please run first example before: 1-simple_test.py")
    exit(1)





print("\nSimka distance exporter allow you to manipulate the results of Simka.")
print("\tBy default, simka provides distance matrices containing all processed datasets, the distance exporter allow you to:")
print("\t1) Quickly create new distance matrices containing only your desired datasets")
print("\t2) Set the order of datasets in the distance matrices")
print("\n\tThe executable " + export_bin + " is used to export the distance matrices, its parameters are the follow:")
print("\t-in: dir to simka distance matrices in binary format. This directory called \"matrix_binary\" is located in simka results dir")
print("\t-out: output directory for new distance matrices")
print("\t-in-ids: a file containing the desired datasets in the new distance matrices (one dataset ID per line)")


#-----------------------------------------------------------------------------
# Export distance matrices with all datasets available (no order specified)
#-----------------------------------------------------------------------------
command = export_bin
command += " -in " + input_dir
command += " -out " + "distance_matrices_all"
os.system(command + " > /dev/null 2>&1  ")

print("\n\nExport distance matrices with all datasets available (no order specified)")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# Export distance matrices with all dataset available (in order A B C D E)
#-----------------------------------------------------------------------------
command = export_bin
command += " -in " + input_dir
command += " -out " + "distance_matrices_ABCDE"
command += " -in-ids " + "data/export_ABCDE.txt"
os.system(command + " > /dev/null 2>&1  ")

print("\nExport distance matrices with all datasets available (in order A B C D E)")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# Export distance matrices with only B and D datasets
#-----------------------------------------------------------------------------
command = export_bin
command += " -in " + input_dir
command += " -out " + "distance_matrices_BD"
command += " -in-ids " + "data/export_BD.txt"
os.system(command + " > /dev/null 2>&1  ")

print("\nExport distance matrices with only B and D datasets")
print("\tcommand used:")
print("\t" + command)
