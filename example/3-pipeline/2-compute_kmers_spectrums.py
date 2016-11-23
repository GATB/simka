
import os, shutil, sys

if len(sys.argv) == 1:
	print("usage: python 2-compute_kmers_spectrums.py <simka2_script_dir> <database_dir> <tmp_dir> <input_filename>")
	print("example: python ../example/2-pipeline/2-compute_kmers_spectrums.py ../scripts/simka2/ ./tmp ./tmp ../example/1-simple_example/simka_input.txt")
	exit(1)

simka2_scripts = sys.argv[1]
database_dir = sys.argv[2]
tmp_dir = sys.argv[3]
input_filename = sys.argv[4]

databaseDir = os.path.join(database_dir, "simka_database")

#-----------------------------------------------------------------------------
# Create dir for temporary computation files
#-----------------------------------------------------------------------------
tmpComputationDir = os.path.join(tmp_dir, "simka_tmp")
if os.path.exists(tmpComputationDir):
	shutil.rmtree(tmpComputationDir)
os.makedirs(tmpComputationDir)


#-----------------------------------------------------------------------------
# Compute k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka2_scripts, "core", "simka2-count.py")
command += " -database-dir " + databaseDir
command += " -out-tmp " + tmpComputationDir
command += " -in " + input_filename

print("Computing k-mer spectrums:")
print("\t" + command)

os.system(command)

