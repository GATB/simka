
import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])



databaseDir = "./simka_database"

simka2_script_dir = "../../scripts/simka2"


#-----------------------------------------------------------------------------
# Compute distances between k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka2_script_dir, "core", "simka2-distance.py")
command += " -database-dir " + databaseDir
os.system(command + " > /dev/null 2>&1  ")

print("\nStep 3) Compute distance between k-mer spectrums")
print("\tcommand used:")
print("\t" + command)

print("\nThe script " + os.path.join(simka2_script_dir, "core", "simka2-distance.py") + " is used to compute the distances")
print("\tThis script must be linked to a simka database (-database-dir)")

print("\nThe distance binaries are store in " + os.path.join(databaseDir, "distance", "matrix_binary"))
print("The distance exporter can be run on those binaries to extract distance matrices (see example 1-basic_usage/3-export-distance)")
#-----------------------------------------------------------------------------
# Run the distance exporter
# The distance exporter reads matrices in binary format (-in)
# extracts the rows/columns depending on supplied datasets id (-in-ids)
# and outputs distance matrices in ascii format readable by R (-out)
#-----------------------------------------------------------------------------
#reprise ici!!!!
#command = os.path.join(SCRIPT_DIR, "bin", "simka2-export") + \
#	" -out " + args.output_dir + \
#		" -in " + matrixBinaryFilename
#os.system(command)