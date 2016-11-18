
import os, sys, argparse, shutil
from simka2_database import SimkaDatabase
from simka2_utils import JobScheduler, Simka2ResourceAllocator, ProgressBar

#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-in', action="store", dest="_inputFilename")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
parser.add_argument('-database-dir', action="store", dest="_databaseDir")
#parser.add_argument('-simka-scripts', action="store", dest="_simkaScriptsDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")

parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

#----------------------------------------------------------------------
#----------------------------------------------------------------------

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]


#----------------------------------------------------------------------
#----------------------------------------------------------------------
class ComputeKmerSpectrumAll():

	def __init__(self):
		self.database = SimkaDatabase(args._databaseDir)
		self.nbDatasetToProcess = 0
		self.processed_ids = []
		#self.jobCores = None
		#self.jobMemory = None

	def execute(self):

		self.countNbDatasetToProcess()

		self.resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs))
		#self.resourceAllocator.maxJobMerge = self.database._nbPartitions
		#self.resourceAllocator.nbSamples = self.nbDatasetToProcess
		maxJobs, self.jobCores, self.jobMemory = self.resourceAllocator.executeForCountJobs(self.nbDatasetToProcess)

		self.jobScheduler = JobScheduler(maxJobs, ProgressBar("Computing k-mer spectrums", self.nbDatasetToProcess))

		self.jobScheduler.start()
		self.computeKmerSpectrums()
		self.jobScheduler.join()

		#self.database.save()

	def computeKmerSpectrums(self):

		inputFile = open(args._inputFilename, "r")
		for line in inputFile:
			line = line.strip()
			if line == "": continue

			id, filenames = line.replace(" ", "").split(":")

			#-- checkpoint
			isDatasetProcessed = self.database.contains_entry(id)

			if isDatasetProcessed:
				continue
			else:
				pass

			#-- create dir for temporary files
			outputDirTemp = self.getOutputDirTemp(id)
			if os.path.exists(outputDirTemp):
				shutil.rmtree(outputDirTemp)

			os.makedirs(outputDirTemp)

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

			self.computeKmerSpectrum(id, inputFilename, outputDirTemp, nbPairedDatasets)



	def computeKmerSpectrum(self, id, inputFilename, outputDirTemp, nbPairedDatasets):

		#print id
		kmerSpectrumOutputDir = os.path.join(self.database.get_default_kmer_spectrum_dir_of_id(id, True))

		#kmerSpectrumDir = os.path.join(self.database.dirname, args._relativeKmerSpectrumDir)

		if os.path.exists(kmerSpectrumOutputDir):
			shutil.rmtree(kmerSpectrumOutputDir)

		os.makedirs(kmerSpectrumOutputDir)

		#if not os.path.exists(args._outputDirTemp):
		#	os.makedirs(args._outputDirTemp)

		#if os.path.exists(outputDir):
		#	shutil.rmtree(outputDir)
		#	os.mkdir(outputDir)

		#if not os.path.exists(args._outputDirTemp):
		#	os.makedirs(args._outputDirTemp)

		#command = "python " + os.path.join(SCRIPT_DIR, "compute_kmer_spectrum.py") + " " + \
		#		  " -in " + inputFilename + \
		#		  " -id " + str(id) + \
		#		  " -out-tmp " + outputDirTemp + \
		#		  " -database " + args._databaseDir + \
		#		  " -out " + outputDir + \
		#		  " -nb-dataset " + str(nbPairedDatasets) + \
		#		  " -simka-bin " + args._simkaBinDir
		command = os.path.join(args._simkaBinDir, "simkaCountProcess") + " " + \
			os.path.join(args._simkaBinDir,"simka2-count") + \
			" -id " + id + \
			" -in " + inputFilename + \
			" -out-tmp " + outputDirTemp + \
			" -out " + kmerSpectrumOutputDir + \
			" -kmer-size " + str(self.database._kmerSize) + \
			" -nb-dataset " + str(nbPairedDatasets) + \
			" -abundance-min " + str(self.database._abundanceMin) + \
			" -nb-partitions " + str(self.database._nbPartitions) + \
			" -max-memory " + str(self.jobMemory) + \
			" -nb-cores " + str(self.jobCores) + \
			"   > /dev/null 2>&1     &"

		#print("compute_kmer_spectrums_all.py: Add log file system")

		print(command)
		os.system(command)

		checkPointFilename = os.path.join(kmerSpectrumOutputDir, "success")
		self.jobScheduler.submitJob((checkPointFilename, self.jobEnd, (id, self.database.get_default_kmer_spectrum_dir_of_id(id, False), outputDirTemp)))

	def jobEnd(self, data):
		id = data[0]
		kmerSpectrumOutputDir = data[1]
		outputDirTemp = data[2]

		self.database.add_entry(id, kmerSpectrumOutputDir)

		shutil.rmtree(outputDirTemp)

		#print "\n", id, "\n"

	def countNbDatasetToProcess(self):

		inputFile = open(args._inputFilename, "r")

		for line in inputFile:
			line = line.strip()
			if line == "": continue

			id, filenames = line.replace(" ", "").split(":")

			#-- checkpoint
			isDatasetProcessed = self.database.contains_entry(id)

			if isDatasetProcessed:
				continue
			else:
				pass

			self.nbDatasetToProcess += 1

		inputFile.close()

	def getOutputDirTemp(self, id):
		return os.path.join(args._outputDirTemp, id + "_temp")



c = ComputeKmerSpectrumAll()
c.execute()