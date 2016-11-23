

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

print("\ncommand used:")
print("\t" + command)
