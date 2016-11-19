
import os, sys, argparse, shutil
from simka2_database import SimkaDatabase
from simka2_utils import Simka2ResourceAllocator, JobScheduler, ProgressBar


#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")

parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")

args =  parser.parse_args()


class SimkaKmerSpectrumMerger():

	#MAX_OPEN_FILES = 1000
	#MAX_OPEN_FILES_PER_MERGE = 100
	MAX_OPEN_FILES = 1000
	MAX_OPEN_FILES_PER_MERGE = 2

	def __init__(self):
		self.database = SimkaDatabase(args._databaseDir)

	def clearTempDir(self):
		self.tempDir = os.path.join(self.database.dirname, "merge")
		if os.path.exists(self.tempDir):
			shutil.rmtree(self.tempDir)
		os.makedirs(self.tempDir)

		print "faire une routine qui check si il y a des merge dir inutile (recup merge id dans database 2) verifier s'il y a des dir non referencer)"

	def execute(self):

		maxJobsByOpenFile = SimkaKmerSpectrumMerger.MAX_OPEN_FILES / SimkaKmerSpectrumMerger.MAX_OPEN_FILES_PER_MERGE
		self.resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs))
		maxJobs, jobCores = self.resourceAllocator.executeForDistanceJobs(self.database._nbPartitions)
		maxJobs = min(maxJobs, maxJobsByOpenFile)

		self.jobScheduler = JobScheduler(maxJobs, ProgressBar("Merging k-mer spectrums", self.database._nbPartitions))

		#---
		self.mergeKmerSpectrums()


	def getNextMergeID(self):
		dirs = [f for f in os.listdir(self.database.mergeKmerSpectrumAbsDir) if os.path.isdir(os.path.join(self.database.mergeKmerSpectrumAbsDir, f))]
		#print("simka2-merge:  ", dirs)
		dirIndex = set()
		for dir in dirs:
			dirIndex.add(int(dir))
		id = 0
		while id in dirIndex:
			id += 1
		return str(id)


	def mergeKmerSpectrums(self):

		"""
		usedDatasetIds = {}

		spectrums = []
		for id in self.database.entries:

			kmerSpectrumDir = self.database.entries_infos[id]
			id = self.database.get_id_from_dir(kmerSpectrumDir)

			if id in usedDatasetIds:
				continue
			usedDatasetIds[id] = True

			spectrums.append(id)


		while(len(spectrums) > SimkaKmerSpectrumMerger.MAX_OPEN_FILES_PER_MERGE):

			mergeDatasetDirs = []

			sys.stdout.write("Merging ")

			for i in range(0, SimkaKmerSpectrumMerger.MAX_OPEN_FILES_PER_MERGE):
				id = spectrums[i]
				mergeDatasetDirs.append(id)
				sys.stdout.write(mergeDatasetDirs[i] + " ")

			mergeDestID = mergeDatasetDirs[0]
			sys.stdout.write( "    in    " + mergeDestID)
			print("")
			#cout << "    in    " << mergeDestID << endl;


			for i in range(0, len(mergeDatasetDirs)):
				del spectrums[0]

			self.clearTempDir()
			self.jobScheduler.start()
			self.mergeSomeKmerSpectrums(os.path.join(self.tempDir, "__merge_input.txt"), mergeDatasetDirs, mergeDestID)
			self.jobScheduler.join()
			self.mergeEnd(os.path.join(self.tempDir, "__merge_input.txt"), mergeDatasetDirs, mergeDestID)

			spectrums.insert(0, mergeDestID)

		"""
		list_kmerSpectrumDirs_Sizes = []
		usedDirs = set()

		for id in self.database.entries:

			kmerSpectrumDir = self.database.get_kmer_spectrum_dir_of_id(id, False)
			#id = self.database.get_id_from_dir(kmerSpectrumDir)

			#id = self.database.get_id_from_dir(kmerSpectrumDir)
			#datasetID = DATASET_IDS[i]
			#print(datasetID + " " + getLinkedDataset(datasetID));

			#datasetID = getLinkedDataset(datasetID)

			if kmerSpectrumDir in usedDirs:
				continue
			usedDirs.add(kmerSpectrumDir)

			kmerSpectrumAbsDir = self.database.get_kmer_spectrum_dir_of_id(id, True)
			size = self.getDirSize(kmerSpectrumAbsDir)

			list_kmerSpectrumDirs_Sizes.append((kmerSpectrumDir, size))

			#print("used: " + id)



		while(len(list_kmerSpectrumDirs_Sizes) > SimkaKmerSpectrumMerger.MAX_OPEN_FILES_PER_MERGE):

			#//cout << "Start merging pass" << endl;
			#sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);
			list_kmerSpectrumDirs_Sizes.sort(key=lambda x: x[1])

			#for i in range(0, len(list_kmerSpectrumDirs_Sizes)):
			#	print i, list_kmerSpectrumDirs_Sizes[i]

			#//vector<size_t> mergeDatasetIds;
			mergeDatasetDirs = []
			#//vector<size_t> toRemoveItem;

			sys.stdout.write("Merging ")

			#//cout << endl;
			#//cout << "merging" << endl;
			for i in range(0, SimkaKmerSpectrumMerger.MAX_OPEN_FILES_PER_MERGE):
				sfi = list_kmerSpectrumDirs_Sizes[i]
				#//mergeDatasetIds.push_back(get<2>(sfi));
				mergeDatasetDirs.append(sfi[0])
				#//cout << get<1>(sfi) << endl;
				#//datasetIndex += 1;
				#//if(datasetIndex >= _nbBanks) break;

				sys.stdout.write(mergeDatasetDirs[i] + " ")
				#cout << getDatasetID(mergeDatasetDirs[i]) << " ";
				#//cout << "First val must never be greater than second:   " << i << "  " << _nbBanks << endl;
				#//cout << "\t" << get<1>(sfi) << endl;
			#}

			#mergeDestID = mergeDatasetDirs[0]
			mergeDestID = self.getNextMergeID()

			sys.stdout.write( "    in    " + mergeDestID)
			print("")
			#cout << "    in    " << mergeDestID << endl;


			for i in range(0, len(mergeDatasetDirs)):
				#_links[getDatasetID(mergeDatasetDirs[i])] = mergeDestID;
				#updateLinks(getDatasetID(mergeDatasetDirs[i]), mergeDestID);

				#filenameSizes.erase(filenameSizes.begin());
				del list_kmerSpectrumDirs_Sizes[0]

			mergeOutputRelativeDir = os.path.join(self.database.mergeKmerSpectrumRelativeDir, mergeDestID)
			mergeOutputAbsDir = os.path.join(self.database.dirname, mergeOutputRelativeDir)
			#mergeOutputDir = os.path.join(self.database.mergeKmerSpectrumDir, mergeDestID)
			os.mkdir(mergeOutputAbsDir)

			self.clearTempDir()
			self.jobScheduler.start()
			self.mergeSmallestKmerSpectrums(os.path.join(self.tempDir, "__merge_input.txt"), mergeDatasetDirs, mergeOutputRelativeDir)
			self.jobScheduler.join()
			self.mergeEnd(os.path.join(self.tempDir, "__merge_input.txt"), mergeDatasetDirs, mergeOutputRelativeDir)

			#print("a changeer, path complet au lieu de juste id dans list")
			list_kmerSpectrumDirs_Sizes.append((mergeOutputRelativeDir, self.getDirSize(mergeOutputAbsDir)))
			#exit(1)


	def mergeSmallestKmerSpectrums(self, merge_input_filename, mergedDirs, mergeOutputRelativeDir):



		#if not os.path.exists()
		#print dataset_to_merge_ids
		#merge_input_filename = os.path.join(self.database.dirname, "__merge_input.txt")
		merge_input_file = open(merge_input_filename, "w")
		for relDir in mergedDirs:
			merge_input_file.write(os.path.join(self.database.dirname, relDir) + "\n")
		merge_input_file.close()

		#print dataset_to_merge_ids

		#print id
		#outputDir = os.path.join(DATABASE.get_kmer_spectrum_dir_of_id(id, True))

		#successFilename = os.path.join(outputDir, "success")

		#if os.path.exists(successFilename):
		#	return

		#if os.path.exists(outputDir):
		#	shutil.rmtree(outputDir)
		#	os.mkdir(outputDir)

		#self.mergeEndData = []
		#self.mergeEndData[0] = merge_input_filename
		#self.mergeEndData[1] = dataset_to_merge_ids
		#self.mergeEndData[2] = merge_dest_id

		#command = "./bin/simkaCountProcess \"./bin/simka2-computeKmerSpectrum -id F1 -in /Users/gbenoit/workspace/gits/gatb-simka/example/B.fasta -out-tmp /local/output/o_simka2_test/ -out /local/output/o_simka2_result\""
		for i in range(0, self.database._nbPartitions):
			command = os.path.join(args._simkaBinDir, "simka2-merge") + \
				" -in " + merge_input_filename + \
				" -database-dir " + args._databaseDir + \
				" -kmer-size " + str(self.database._kmerSize) + \
				" -partition-id " + str(i) + \
				" -out " + os.path.join(self.database.dirname, mergeOutputRelativeDir) + \
				"   > /dev/null 2>&1     &"
			print command
			os.system(command)

			checkPointFilename = os.path.join(self.tempDir, str(i) + "-success")
			self.jobScheduler.submitJob((checkPointFilename, self.jobEnd, ()))


	def jobEnd(self, data):
		pass

	def mergeEnd(self, merge_input_filename, mergedDirs, mergeOutputRelativeDir):

		#mergeDir = self.database.get_kmer_spectrum_dir_of_id(merge_dest_id, True)

		#for i in range(0, self.database._nbPartitions):
		#	os.remove(os.path.join(mergeDir, str(i) + ".gz"))
		#	shutil.move(os.path.join(mergeDir, str(i) + ".gz.temp"), os.path.join(mergeDir, str(i) + ".gz"))


		#save merge infos
		command = os.path.join(args._simkaBinDir, "simka2-merge") + \
			" -in " + merge_input_filename + \
			" -database-dir " + args._databaseDir + \
			" -kmer-size " + str(self.database._kmerSize) + \
			" -partition-id " + "0" + \
			" -out " + os.path.join(self.database.dirname, mergeOutputRelativeDir) + \
			" -save-merge-info " + \
			"   > /dev/null 2>&1"
		ret = os.system(command)
		if ret != 0: exit(1)

		dirs_to_delete = []
		for relDir in mergedDirs:
			dir_that_have_been_merged = os.path.join(self.database.dirname, relDir)
			dirs_to_delete.append(dir_that_have_been_merged)

		#dataset_to_merge_ids.remove(merge_dest_id)

		dataset_to_merge_ids_index = set()
		for relDir in mergedDirs:
			dataset_to_merge_ids_index.add(relDir)
		self.database.change_entries(dataset_to_merge_ids_index, mergeOutputRelativeDir)
		self.database.save()


		for dir in dirs_to_delete:
			shutil.rmtree(dir)

		#exit(1)



	def getDirSize(self, dirpath):
		return sum(os.path.getsize(os.path.join(dirpath, f)) for f in os.listdir(dirpath))


c = SimkaKmerSpectrumMerger()
c.execute()

