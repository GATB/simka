
import os, sys, argparse, shutil
from simka_database import SimkaDatabase



#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
#parser.add_argument('-in', action="store", dest="_inputFilename")
#parser.add_argument('-out', action="store", dest="_outputDir")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

#----------------------------------------------------------------------
#----------------------------------------------------------------------

DATABASE = SimkaDatabase(args._databaseDir)
SIMKA2_SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

#NB_BANKS_TO_MERGE = 0
#DATASET_IDS = []
#----------------------------------------------------------------------
#----------------------------------------------------------------------

def mergeKmerSpectrums():
	merge_script_filename = os.path.join(SIMKA2_SCRIPT_DIR, "./simka2-merge.py")
	command = "python " + merge_script_filename + " -database-dir " + DATABASE.dirname + " -simka-bin " + args._simkaBinDir
	os.system(command)

def computeDistance():

	#uniqKmerSpectrumDirs = set()

	#for id, kmerSpectrumDir in DATABASE.entries_infos.items():
	#	uniqKmerSpectrumDirs.add(kmerSpectrumDir)

	#distanceInputFilename = os.path.join(DATABASE.dirname, "distance", "_input.txt")
	#distanceInputFile = open(distanceInputFilename, "w")
	#for dir in uniqKmerSpectrumDirs:
	#	distanceInputFile.write(dir + "\n")
	#distanceInputFile.close()

	for i in range(0, DATABASE._nbPartitions):
		command = os.path.join(args._simkaBinDir, "simka2-distance") + \
			" -database-dir " + args._databaseDir + \
			" -kmer-size " + str(DATABASE._kmerSize) + \
			" -partition-id " + str(i)
		print command
		os.system(command)
		#print("lala")
		#break

	command = os.path.join(args._simkaBinDir, "simka2-distance-final") + \
		" -database-dir " + args._databaseDir + \
		" -kmer-size " + str(DATABASE._kmerSize) + \
		" -nb-partitions " + str(DATABASE._nbPartitions)
	os.system(command)


#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():

	print ("Attention remettere le merge dans simka2-distance main()")
	#mergeKmerSpectrums()
	computeDistance()



main()
