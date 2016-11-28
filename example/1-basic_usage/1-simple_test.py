

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])





#-----------------------------------------------------------------------------
# Set required simka parameters
#-----------------------------------------------------------------------------
simka_script_dir = "../../scripts/simka2"
input_filename = "../data/simka_input.txt"
output_dir = "simka_results"
temp_dir = "simka_temp"


#-----------------------------------------------------------------------------
# Create the command to execute
#-----------------------------------------------------------------------------
command = "python "
command += simka_script_dir + "/simka.py"
command += " -in " + input_filename
command += " -out " + output_dir
command += " -out-tmp " + temp_dir

#-----------------------------------------------------------------------------
# Run Simka
#-----------------------------------------------------------------------------
os.system(command + " > /dev/null 2>&1  ")

#-----------------------------------------------------------------------------
# Clear temporary computation dir
#-----------------------------------------------------------------------------
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

print("\nthis is a simple example of running simka, required parameter are:")
print("\t-in: file containing the paths of dataset to process and their ID (one dataset per line)")
print("\t\tID1: filename.fasta")
print("\t\tID2: filename2.fasta")
print("\t\t...")
print("\t-out-tmp: dir for simka temporary computation files. Target your faster available disk")
print("\t-out (optional): output directory for simka results")
print("\ncommand used:")
print("\t" + command)
