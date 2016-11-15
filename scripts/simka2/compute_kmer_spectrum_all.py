
import os, sys, argparse, shutil
from simka_database import SimkaDatabase


#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-in', action="store", dest="_inputFilename")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
parser.add_argument('-database', action="store", dest="_databaseDir")
#parser.add_argument('-simka-scripts', action="store", dest="_simkaScriptsDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

#----------------------------------------------------------------------
#----------------------------------------------------------------------

DATABASE = SimkaDatabase(args._databaseDir)
SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

if not os.path.exists(args._outputDirTemp):
	os.makedirs(args._outputDirTemp)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

def getOutputDirTemp(id):
	return os.path.join(args._outputDirTemp, id + "_temp")

def computeKmerSpectrum(id, inputFilename, outputDirTemp, nbPairedDatasets):

	#print id
	outputDir = os.path.join(DATABASE.get_default_kmer_spectrum_dir_of_id(id, True))

	#if os.path.exists(outputDir):
	#	shutil.rmtree(outputDir)
	#	os.mkdir(outputDir)

	command = "python " + os.path.join(SCRIPT_DIR, "compute_kmer_spectrum.py") + " " + \
			  " -in " + inputFilename + \
			  " -id " + str(id) + \
			  " -out-tmp " + outputDirTemp + \
			  " -database " + args._databaseDir + \
			  " -out " + outputDir + \
			  " -nb-dataset " + str(nbPairedDatasets) + \
			  " -simka-bin " + args._simkaBinDir

	#command = "./bin/simkaCountProcess \"./bin/simka2-computeKmerSpectrum -id F1 -in /Users/gbenoit/workspace/gits/gatb-simka/example/B.fasta -out-tmp /local/output/o_simka2_test/ -out /local/output/o_simka2_result\""
	#command = SIMKA_BIN + "simka2-computeKmerSpectrum -id " + id + " -in " + filenames + " -out-tmp " + args._outputDirTemp + " -out " + outputDir + " -kmer-size " + str(DATABASE._kmerSize)
	print(command)
	os.system(command)
	DATABASE.add_entry(id, os.path.join(DATABASE.get_default_kmer_spectrum_dir_of_id(id, False)))

#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
	inputFile = open(args._inputFilename, "r")
	for line in inputFile:
		line = line.strip()
		if line == "": continue

		id, filenames = line.replace(" ", "").split(":")

		#-- checkpoint
		isDatasetProcessed = DATABASE.contains_entry(id)

		#print id, isDatasetProcessed
		if isDatasetProcessed:
			continue
			#outputDir = os.path.join(DATABASE.get_kmer_spectrum_dir_of_id(id, True))
			#successFilename = os.path.join(outputDir, "success")
			#if os.path.exists(successFilename):
			#	return
		else:
			pass

		#-- create dir for temporary files
		outputDirTemp = getOutputDirTemp(id)
		if os.path.exists(outputDirTemp):
			shutil.rmtree(outputDirTemp)

		os.mkdir(outputDirTemp)

		#-- create input filename
		inputFilename = os.path.join(outputDirTemp, "input.txt")
		inputFile = open(inputFilename, "w")

		abs_filenames = ""
		filenames_paireds = filenames.split(";")
		nbPairedDatasets = len(filenames_paireds)
		for filenames_paired in filenames_paireds:
			filenames_paireds_concats = filenames_paired.split(",")
			for filenames_paireds_concat in filenames_paireds_concats:
				if filenames_paireds_concat[0] == "/":
					abs_filenames +=  filenames_paireds_concat
				else:
					dir = os.path.dirname(os.path.realpath(args._inputFilename))
					abs_filenames +=  dir + "/" + filenames_paireds_concat
				abs_filenames += "\n"
			abs_filenames = abs_filenames[:-1]
			abs_filenames += "\n"
		abs_filenames = abs_filenames[:-1]

		inputFile.write(abs_filenames)
		inputFile.close()

		computeKmerSpectrum(id, inputFilename, outputDirTemp, nbPairedDatasets)

main()
DATABASE.save()