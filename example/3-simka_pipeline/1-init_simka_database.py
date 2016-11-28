
import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])




if len(sys.argv) == 1:
	print("usage: python 1-init_simka_database.py <kmer_size>")
	print("example: python ../example/3-simka_pipeline/1-init_simka_database.py 31")
	exit(1)

kmer_size = sys.argv[1]
databaseDir = "./simka_database"


simka2_script_dir = "../../scripts/simka2"

#-----------------------------------------------------------------------------
# Create database dir
#-----------------------------------------------------------------------------
if os.path.exists(databaseDir):
	shutil.rmtree(databaseDir)
os.makedirs(databaseDir)

#-----------------------------------------------------------------------------
# Init Simka database
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka2_script_dir, "core", "simka2-init.py")
command += " -database-dir " + databaseDir
command += " -kmer-size " + kmer_size
os.system(command + " > /dev/null 2>&1  ")

print("\nStep 1) Init simka database")
print("\tcommand used:")
print("\t" + command)

print("\nYou have defined the parameters that Simka will used to compute k-mer spectrums and distances")
