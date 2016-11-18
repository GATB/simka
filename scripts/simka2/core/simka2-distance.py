
import os, sys, argparse, shutil
from simka2_database import SimkaDatabase
from simka2_utils import Simka2ResourceAllocator, JobScheduler, ProgressBar

parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
#parser.add_argument('-in', action="store", dest="_inputFilename")
#parser.add_argument('-out', action="store", dest="_outputDir")

parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

SIMKA2_SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

#----------------------------------------------------------------------
#----------------------------------------------------------------------

class Simka_ComputeDistance():

	def __init__(self):
		self.database = SimkaDatabase(args._databaseDir)
		self.clearTempDir()

	def clearTempDir(self):
		self.tempDir = os.path.join(self.database.dirname, "distance", "temp_parts")
		if os.path.exists(self.tempDir):
			shutil.rmtree(self.tempDir)
		os.makedirs(self.tempDir)

	def execute(self):


		#print ("Attention remettere le merge dans simka2-distance main()")
		self.mergeKmerSpectrums()

		self.resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs))
		maxJobs, jobCores = self.resourceAllocator.executeForDistanceJobs(self.database._nbPartitions)

		self.jobScheduler = JobScheduler(maxJobs, ProgressBar("Computing distances", self.database._nbPartitions))

		self.jobScheduler.start()
		self.computeDistanceParts()
		self.jobScheduler.join()

		self.computeDistanceFinal()

	def mergeKmerSpectrums(self):
		merge_script_filename = os.path.join(SIMKA2_SCRIPT_DIR, "./simka2-merge.py")
		command = "python " + merge_script_filename + " -database-dir " + self.database.dirname + " -simka-bin " + args._simkaBinDir
		os.system(command)

	def computeDistanceParts(self):

		#uniqKmerSpectrumDirs = set()

		#for id, kmerSpectrumDir in DATABASE.entries_infos.items():
		#	uniqKmerSpectrumDirs.add(kmerSpectrumDir)

		#distanceInputFilename = os.path.join(DATABASE.dirname, "distance", "_input.txt")
		#distanceInputFile = open(distanceInputFilename, "w")
		#for dir in uniqKmerSpectrumDirs:
		#	distanceInputFile.write(dir + "\n")
		#distanceInputFile.close()

		for i in range(0, self.database._nbPartitions):
			self.computeDistancePart(i)

	def computeDistancePart(self, partitionId):
		command = os.path.join(args._simkaBinDir, "simka2-distance") + \
			" -database-dir " + args._databaseDir + \
			" -kmer-size " + str(self.database._kmerSize) + \
			" -partition-id " + str(partitionId) + \
			"   > /dev/null 2>&1     &"
		print command
		os.system(command)

		checkPointFilename = os.path.join(self.tempDir, str(partitionId) + "-success")
		self.jobScheduler.submitJob((checkPointFilename, self.jobEnd, ()))
		#print("lala")
		#break

	def computeDistanceFinal(self):
		command = os.path.join(args._simkaBinDir, "simka2-distanceFinal") + \
			" -database-dir " + args._databaseDir + \
			" -kmer-size " + str(self.database._kmerSize) + \
			" -nb-partitions " + str(self.database._nbPartitions)
		os.system(command)

	def jobEnd(self, data):
		pass
		#id = data[0]
		#kmerSpectrumOutputDir = data[1]
		#outputDirTemp = data[2]

		#self.database.add_entry(id, kmerSpectrumOutputDir)

		#shutil.rmtree(outputDirTemp)

		#print "\n", id, "\n"


c = Simka_ComputeDistance()
c.execute()
