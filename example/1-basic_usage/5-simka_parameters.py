

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])


#-----------------------------------------------------------------------------
# Set required simka parameters
#-----------------------------------------------------------------------------
simka_script_dir = "../../scripts/simka2"
input_filename = "../data/simka_input.txt"
temp_dir = "simka_temp"

def initSimkaCommand():
    command = "python "
    command += simka_script_dir + "/simka.py"
    command += " -in " + input_filename
    command += " -out-tmp " + temp_dir
    return command

#-----------------------------------------------------------------------------
# Core parameters
#-----------------------------------------------------------------------------
command = initSimkaCommand()
command += " -nb-cores " + "4"
command += " -max-memory " + "4000"
#os.system(command + " > /dev/null 2>&1  ")

print("\n\nCore options: ")
print("\tCores options increase the performance of Simka.")
print("\t-nb-cores X: Simka will use X cores to speed up computation")
print("\t-max-memory X: Simka will use up to X MB of memory during computation")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# Distance parameters
#-----------------------------------------------------------------------------
command = initSimkaCommand()
command += " -simple-dist "
command += " -complex-dist "
#os.system(command + " > /dev/null 2>&1  ")

print("\n\nDistance options: ")
print("\tUse distance options if you want more distance information. But be aware that Simka will need more time (and memory) to compute them.")
print("\tBy default, Simka provide abundance-based Bray-Curtis distance and Presence-Absence-based Jaccard distance")
print("\t-simple-dist: simka will compute several distance that are fast to compute")
print("\t-complex-dist: simka will compute several distances that can be very long to compute")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# K-mer parameters
#-----------------------------------------------------------------------------
command = initSimkaCommand()
command += " -kmer-size " + "31"
command += " -abundance-min " + "2"
command += " -abundance-max " + "10000"
#os.system(command + " > /dev/null 2>&1  ")

print("\n\nK-mer options: ")
print("\tUse k-mer options to indicate the size of k-mer considered and filter based on k-mer\'s abundance")
print("\t-kmer-size X: Simka will use k-mer of size X")
print("\t-abundance-min X: Simka will discard all distinct k-mers with count < X in a given dataset")
print("\t\tBy default, Simka discards all k-mer seen only one time (-abundance-min 2).")
print("\t\tIf you have dataset with low coverage it can be good to disable this filter by setting -abundance-min to 0 or 1")
print("\t\tBut be aware that the time and ressources required by Simka can be a lot more if you disable it.")
print("\t-abundance-max X: Simka will discard all distinct k-mers with count > X in a given dataset")

print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# K-mer parameters
#-----------------------------------------------------------------------------
command = initSimkaCommand()
command += " -max-reads " + "100000"
command += " -min-read-size " + "70"
command += " -min-shannon-index " + "1.5"
#os.system(command + " > /dev/null 2>&1  ")

print("\n\nRead options: ")
print("\tRead options are quality control options at the read level")
print("\t-max-reads X: simka will use the X first reads in each dataset")
print("\t-min-read-size X: Simka will discard any read with size < X")
print("\t-min-shannon-index X: Simka will discard any read with a sequence complexity < X")

print("\tcommand used:")
print("\t" + command)