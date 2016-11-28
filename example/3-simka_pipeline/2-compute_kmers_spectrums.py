
import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])




input_filename = "../data/simka_input.txt"
databaseDir = "./simka_database"
tmpComputationDir = "./simka_temp"

if not os.path.exists(databaseDir):
	print("Please run first example before: 1-init_simka_database.py")
	exit(1)

simka2_script_dir = "../../scripts/simka2"

#-----------------------------------------------------------------------------
# Create dir for temporary computation files
#-----------------------------------------------------------------------------
if os.path.exists(tmpComputationDir):
	shutil.rmtree(tmpComputationDir)
os.makedirs(tmpComputationDir)


#-----------------------------------------------------------------------------
# Compute k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka2_script_dir, "core", "simka2-count.py")
command += " -database-dir " + databaseDir
command += " -out-tmp " + tmpComputationDir
command += " -in " + input_filename
os.system(command + " > /dev/null 2>&1  ")

print("\nStep 2) Compute k-mer spectrums")
print("\tcommand used:")
print("\t" + command)

print("\nYou have computed the k-mer spectrum of each input dataset")
print("The k-mer spectrum is the set of each distinct k-mer associated with their abundance (number of ocurences in the dataset)")

print("\nThe script " + os.path.join(simka2_script_dir, "core", "simka2-count.py") + " is used to compute those spectrums")
print("\tThis script must be linked to a simka database (-database-dir)")
print("\tA directory must be specified for temporary computation files (-out-tmp)")
print("\tYou indicate the input file of datasets with the -in option")
