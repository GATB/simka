
import os, shutil, sys

if len(sys.argv) == 1:
	print("usage: python 1-init_simka_database.py <simka2_script_dir> <database_dir> <kmer_size>")
	print("example: python ../example/2-pipeline/1-init_simka_database.py ../scripts/simka2/ ./tmp 31")
	exit(1)

simka2_scripts = sys.argv[1]
database_dir = sys.argv[2]
kmer_size = sys.argv[3]

databaseDir = os.path.join(database_dir, "simka_database")

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
command += os.path.join(simka2_scripts, "core", "simka2-init.py")
command += " -database-dir " + databaseDir
command += " -kmer-size " + kmer_size

print("Creating simka database:")
print("\t" + command)

os.system(command)