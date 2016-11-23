
import os, sys, argparse, shutil
from simka2_database import SimkaDatabase
from simka2_utils import Simka2ResourceAllocator, JobScheduler, ProgressBar, SimkaCommand

parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
#parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTemp")
#parser.add_argument('-in', action="store", dest="_inputFilename")
#parser.add_argument('-out', action="store", dest="_outputDir")

parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")
parser.add_argument('-submit-command', action="store", dest="submit_command", help="command used to submit job")
parser.add_argument('-submit-file', action="store", dest="submit_file", help="filename to a job file template, for HPC system that required a job file")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

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

		self.resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs), args.submit_command, args.submit_file)
		maxJobs, self.jobCores = self.resourceAllocator.executeForDistanceJobs(self.database._nbPartitions)

		self.jobScheduler = JobScheduler(maxJobs, ProgressBar("Computing distances", self.database._nbPartitions))

		self.jobScheduler.start()
		self.computeDistanceParts()
		self.jobScheduler.join()

		self.computeDistanceFinal()

	def mergeKmerSpectrums(self):
		merge_script_filename = os.path.join(SCRIPT_DIR, "./simka2-merge.py")
		command = "python " + merge_script_filename + " -database-dir " + self.database.dirname# + " -simka-bin " + args._simkaBinDir
		command = SimkaCommand.addHPCargs(command, args)
		print command
		ret = os.system(command)
		if ret != 0: exit(1)


	def computeDistanceParts(self):

		#uniqKmerSpectrumDirs = set()

		#for id, kmerSpectrumDir in DATABASE.entries_infos.items():
		#	uniqKmerSpectrumDirs.add(kmerSpectrumDir)

		#distanceInputFilename = os.path.join(DATABASE.dirname, "distance", "_input.txt")
		#distanceInputFile = open(distanceInputFilename, "w")
		#for dir in uniqKmerSpectrumDirs:
		#	distanceInputFile.write(dir + "\n")
		#distanceInputFile.close()
		jobCommandsId = 0

		i = 0
		while(i < self.database._nbPartitions):
			self.jobCommandsFile = open(os.path.join(self.tempDir, "commands_input_" + str(jobCommandsId)), "w")
			for j in range(0, self.jobCores):
				self.constructJobCommand(i)
				i += 1
				if i == self.database._nbPartitions:
					break
			jobCommandsId += 1
			self.jobCommandsFile.close()

		for i in range(0, jobCommandsId):
			checkPointFilename = os.path.join(self.tempDir, "commands_" + str(i) + "-")
			unsuccessCheckPointFilename = checkPointFilename + "unsuccess"
			if os.path.exists(unsuccessCheckPointFilename): os.remove(unsuccessCheckPointFilename)

			command = " python " + os.path.join(SCRIPT_DIR, "simka2-run-job.py") + " "
			command += checkPointFilename + " "
			command += " python " + os.path.join(SCRIPT_DIR, "simka2-run-job-multi.py") + " "
			command += " " + os.path.join(self.tempDir, "commands_input_" + str(i))
			#command += "   > /dev/null 2>&1     &"
			command += " & "
			command = SimkaCommand.createHPCcommand(command, args._isHPC, args.submit_command)
			os.system(command)
			#print command

			self.jobScheduler.submitJob((checkPointFilename, self.jobEnd, ()))
			#print("lala")
			#break


	def constructJobCommand(self, partitionId):

		checkPointFilename = os.path.join(self.tempDir, str(partitionId) + "-")
		unsuccessCheckPointFilename = checkPointFilename + "unsuccess"
		if os.path.exists(unsuccessCheckPointFilename): os.remove(unsuccessCheckPointFilename)

		command = "python " + os.path.join(SCRIPT_DIR, "simka2-run-job.py") + " "
		command += checkPointFilename + " "
		command += os.path.join(SCRIPT_DIR, "..", "bin", "simka2-distance")
		command += " -database-dir " + args._databaseDir
		command += " -kmer-size " + str(self.database._kmerSize)
		command += " -partition-id " + str(partitionId)
		command += "   > /dev/null 2>&1     &"

		self.jobCommandsFile.write(checkPointFilename + "|" + command + "\n")


	def computeDistanceFinal(self):
		command = os.path.join(SCRIPT_DIR, "..", "bin", "simka2-distanceFinal") + \
			" -database-dir " + args._databaseDir + \
			" -kmer-size " + str(self.database._kmerSize) + \
			" -nb-partitions " + str(self.database._nbPartitions)
		print "simka2-distance (computeDistanceFinal):    besoin de pouvoir submit ce job en HPC mode ?"
		#command = SimkaCommand.createHPCcommand(command, args._isHPC, args.submit_command)
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
