
import os, sys, argparse, shutil
from simka_database import SimkaDatabase



#----------------------------------------------------------------------
#----------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir")
#parser.add_argument('-in', action="store", dest="_inputFilename")
#parser.add_argument('-out', action="store", dest="_outputDir")

#print parser.parse_args(['-in', '-bval', '-c', '3'])
args =  parser.parse_args()
#print(args)

#SIMKA_BIN = os.path.dirname(os.path.realpath(__file__)) + "/../../build/bin/"

#----------------------------------------------------------------------
#----------------------------------------------------------------------

DATABASE = SimkaDatabase(args._databaseDir)
#NB_BANKS_TO_MERGE = 0
#DATASET_IDS = []
#----------------------------------------------------------------------
#----------------------------------------------------------------------


def mergeKmerSpectrums(dataset_to_merge_ids, merge_dest_id):

	#print dataset_to_merge_ids
	merge_input_filename = os.path.join(DATABASE.dirname, "__merge_input.txt")
	merge_input_file = open(merge_input_filename, "w")
	for id in dataset_to_merge_ids:
		merge_input_file.write(DATABASE.get_kmer_spectrum_dir_of_id(id, True) + "\n")
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
	
	#command = "./bin/simkaCountProcess \"./bin/simka2-computeKmerSpectrum -id F1 -in /Users/gbenoit/workspace/gits/gatb-simka/example/B.fasta -out-tmp /local/output/o_simka2_test/ -out /local/output/o_simka2_result\""
	for i in range(0, DATABASE._nbPartitions):
		command = os.path.join(args._simkaBinDir, "simka2-merge") + \
				  " -in " + merge_input_filename + \
				  " -kmer-size " + str(DATABASE._kmerSize) + \
				  " -partition-id " + str(i)

		os.system(command)

	#save merge infos
	command = os.path.join(args._simkaBinDir, "simka2-merge") + \
		" -in " + merge_input_filename + \
		" -kmer-size " + str(DATABASE._kmerSize) + \
		" -partition-id " + "0" + \
		" -save-merge-info "
	os.system(command)

	dataset_to_merge_ids.remove(merge_dest_id)
	for id in dataset_to_merge_ids:
		old_dir = DATABASE.get_kmer_spectrum_dir_of_id(id, True)
		shutil.rmtree(old_dir)

		#newRelativeDir = DATABASE.get_kmer_spectrum_dir_of_id(merge_dest_id, False)
		DATABASE.change_entries(id, merge_dest_id)
		#DATABASE.items[id] = newRelativeDir
		#merge_input_file.write( + "\n")


def getDirSize(dirpath):
	return sum(os.path.getsize(os.path.join(dirpath, f)) for f in os.listdir(dirpath))

#def getLinkedDataset(id):
#	return os.path.basename(DATABASE.get_kmer_spectrum_dir_of_id(id, False))




def merge():

	list_kmerSpectrumDirs_Sizes = []
	usedDatasetIds = {}

	for id_dummy, kmerSpectrumDir in DATABASE.entries_infos.items():

		id = DATABASE.get_id_from_dir(kmerSpectrumDir)
		#datasetID = DATASET_IDS[i]
		#print(datasetID + " " + getLinkedDataset(datasetID));

		#datasetID = getLinkedDataset(datasetID)

		if id in usedDatasetIds:
			continue
		usedDatasetIds[id] = True

		kmerSpectrumDir = DATABASE.get_kmer_spectrum_dir_of_id(id, True)
		size = getDirSize(kmerSpectrumDir)

		list_kmerSpectrumDirs_Sizes.append((id, size))

		print("used: " + id)
			#datasetID = _links[datasetID];
			#//cout << datasetID << endl;
			#if (find (usedDatasetIds.begin(), usedDatasetIds.end(), datasetID) != usedDatasetIds.end()){
			#	continue;
			#}
			#usedDatasetIds.push_back(datasetID);

			#string linkedFilename = System::file().getDirectory(_kmerSpectrumDirs[i]) + "/" + datasetID;

			#//filenames.push_back(filename);
			#filenameSizes.push_back(sortItem_Size_Filename_ID(getDatasetSize(linkedFilename), linkedFilename));

			#cout << "used: " << datasetID << endl;
		#}


		#vector<string> partFilenames;


	nbPartLimit = 2

	while(len(list_kmerSpectrumDirs_Sizes) > nbPartLimit):

		#//cout << "Start merging pass" << endl;
		#sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);
		list_kmerSpectrumDirs_Sizes.sort(key=lambda x: x[1])

		#for i in range(0, len(list_kmerSpectrumDirs_Sizes)):
		#	print i, list_kmerSpectrumDirs_Sizes[i]

		#//vector<size_t> mergeDatasetIds;
		mergeDatasetDirs = [];
		#//vector<size_t> toRemoveItem;

		sys.stdout.write("Merging ")

		#//cout << endl;
		#//cout << "merging" << endl;
		for i in range(0, nbPartLimit):
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

		mergeDestID = mergeDatasetDirs[0]
		sys.stdout.write( "    in    " + mergeDestID)
		print("")
		#cout << "    in    " << mergeDestID << endl;


		for i in range(0, len(mergeDatasetDirs)):
			#_links[getDatasetID(mergeDatasetDirs[i])] = mergeDestID;
			#updateLinks(getDatasetID(mergeDatasetDirs[i]), mergeDestID);

			#filenameSizes.erase(filenameSizes.begin());
			del list_kmerSpectrumDirs_Sizes[0]

		mergeKmerSpectrums(mergeDatasetDirs, mergeDestID)

		list_kmerSpectrumDirs_Sizes.append((mergeDestID, getDirSize(DATABASE.get_kmer_spectrum_dir_of_id(mergeDestID, True))))

#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
	pass
	#inputFile = open(args._inputFilename, "r")
	dataset_to_merge_ids = []

	#for line in inputFile:
	#	line = line.strip()
	#	if line == "": continue

	#	id = line.replace(" ", "")
	#	DATASET_IDS.append(id)
		#dataset_to_merge_ids.append(id)

	#NB_BANKS_TO_MERGE = len(dataset_to_merge_ids)
	#mergeKmerSpectrums(dataset_to_merge_ids)



main()
merge()

DATABASE.save()