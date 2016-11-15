
import os, sys, argparse, shutil
from simka_database import SimkaDatabase


#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-in', action="store", dest="_inputFilename")
parser.add_argument('-id', action="store", dest="_datasetID")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
parser.add_argument('-database', action="store", dest="_databaseDir")
parser.add_argument('-out', action="store", dest="_relativeKmerSpectrumDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")
parser.add_argument('-nb-dataset', action="store", dest="_nbPairedDataset")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

#----------------------------------------------------------------------
#----------------------------------------------------------------------

DATABASE = SimkaDatabase(args._databaseDir)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

def computeKmerSpectrum():
	
	id = args._datasetID
	kmerSpectrumDir = os.path.join(DATABASE.dirname, args._relativeKmerSpectrumDir)
	
	if os.path.exists(kmerSpectrumDir):
		shutil.rmtree(kmerSpectrumDir)

	os.makedirs(kmerSpectrumDir)
	
	if not os.path.exists(args._outputDirTemp):
		os.makedirs(args._outputDirTemp)
	
	#print id
	#outputDir = os.path.join(DATABASE.get_kmer_spectrum_dir_of_id(id, True))

#successFilename = os.path.join(outputDir, "success")

#if os.path.exists(successFilename):
#		return

#if os.path.exists(outputDir):
#		shutil.rmtree(outputDir)
#		os.mkdir(outputDir)

	#command = "./bin/simkaCountProcess \"./bin/simka2-computeKmerSpectrum -id F1 -in /Users/gbenoit/workspace/gits/gatb-simka/example/B.fasta -out-tmp /local/output/o_simka2_test/ -out /local/output/o_simka2_result\""
	command = os.path.join(args._simkaBinDir,"simka2-computeKmerSpectrum") + \
		" -id " + id + \
		" -in " + args._inputFilename + \
		" -out-tmp " + args._outputDirTemp + \
		" -out " + kmerSpectrumDir + \
		" -kmer-size " + str(DATABASE._kmerSize) + \
		" -nb-dataset " + args._nbPairedDataset + \
		" -nb-partitions " + str(DATABASE._nbPartitions)

	print(command)
	os.system(command)

	finishFilename = os.path.join(kmerSpectrumDir, "success")
	checkpointFile = open(finishFilename, "w")
	checkpointFile.close()


#DATABASE.add_entry(id, os.path.join(DATABASE.get_kmer_spectrum_dir_of_id(id, False)))

#----------------------------------------------------------------------
#----------------------------------------------------------------------
computeKmerSpectrum()
#DATABASE.save()