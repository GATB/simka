
import os, shutil, sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

simka_script_dir = "../../scripts/simka2"
output_dir = "simka_results"
temp_dir = "simka_temp"

# Clear dirs
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

#-----------------------------------------------------------------------------
# Run simka like in 1-simple_test.py, but keep simka database
#-----------------------------------------------------------------------------
command = "python "
command += simka_script_dir + "/simka.py"
command += " -in " + "../data/simka_input.txt"
command += " -out " + output_dir
command += " -out-tmp " + temp_dir
command += " -keep-tmp " #Keep simka database !
os.system(command + " > /dev/null 2>&1  ")

print("\nRun simka like in 1-simple_test.py, but keep simka database")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# Run simka like in 1-simple_test.py, but with another input file of datasets
#-----------------------------------------------------------------------------
command = "python "
command += simka_script_dir + "/simka.py"
command += " -in " + "../data/simka_input_2.txt"
command += " -out " + output_dir
command += " -out-tmp " + temp_dir
command += " -keep-tmp " #Keep simka database !
os.system(command + " > /dev/null 2>&1  ")

print("\nRun simka like in 1-simple_test.py, but with another input file of datasets (F and G)")
print("\tcommand used:")
print("\t" + command)


print("\nNow, exported matrices in " + output_dir + " includes F and G datasets")