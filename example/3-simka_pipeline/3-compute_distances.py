
import os, shutil, sys

if len(sys.argv) == 1:
	print("usage: python 2-compute_kmers_spectrums.py <simka2_script_dir> <database_dir> <output_dir>")
	print("example: python ../example/2-pipeline/2-compute_kmers_spectrums.py ../scripts/simka2/ ./tmp ./tmp ../example/1-simple_example/simka_input.txt")
	exit(1)

simka2_scripts = sys.argv[1]
database_dir = sys.argv[2]
output_dir = sys.argv[3]

databaseDir = os.path.join(database_dir, "simka_database")

#-----------------------------------------------------------------------------
# Create dir for simka results
#-----------------------------------------------------------------------------
output_dir = os.path.join(output_dir, "simka_results")
if os.path.exists(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir)

#-----------------------------------------------------------------------------
# Compute distances between k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka2_scripts, "core", "simka2-distance.py")
command += " -database-dir " + databaseDir
command += " "
os.system(command)

print("Computing distance between k-mer spectrums:")
print("\t" + command)

os.system(command)

#-----------------------------------------------------------------------------
# Run the distance exporter
# The distance exporter reads matrices in binary format (-in)
# extracts the rows/columns depending on supplied datasets id (-in-ids)
# and outputs distance matrices in ascii format readable by R (-out)
#-----------------------------------------------------------------------------
#reprise ici!!!!
command = os.path.join(SCRIPT_DIR, "bin", "simka2-export") + \
	" -out " + args.output_dir + \
		" -in " + matrixBinaryFilename
os.system(command)