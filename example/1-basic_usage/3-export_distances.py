

import os, shutil, sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

export_bin = "../../scripts/simka2/bin/simka2-export"
input_dir = "simka_results/matrix_binary"

if not os.path.exists(input_dir):
    print("Please run first example before: 1-simple_test.py")
    exit(1)

#-----------------------------------------------------------------------------
# Export distance matrices with all datasets available (no order specified)
#-----------------------------------------------------------------------------
command = export_bin
command += " -in " + input_dir
command += " -out " + "distance_matrices_all"
os.system(command + " > /dev/null 2>&1  ")

print("\nExport distance matrices with all datasets available (no order specified)")
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
