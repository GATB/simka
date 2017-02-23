/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
#include "../utils/SimkaIoUtils.hpp"
#include "../core/SimkaUtils.hpp"
#include "Simka2Utils.hpp"
#include "../minikc/MiniKC.hpp"
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include "Simka2Database.hpp"

//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"



/*
u_int64_t getFileSize(const string& filename){
	std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
	u_int64_t size = in.tellg();
	in.close();
	return size;
}

u_int64_t getDatasetSize(const string& kmerSpectrumDir){
	u_int64_t datasetSize = 0;
	for(size_t partitionID=0; partitionID<SIMKA2_NB_PARTITIONS; partitionID++){
		string partFilename = kmerSpectrumDir + "/" + Stringify::format("%i", partitionID) + ".gz";
		datasetSize += getFileSize(partFilename);
	}
	return datasetSize;
}*/










template<size_t span>
class DatasetMergerWriter : public DiskBasedMergeSort<span>
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;

    struct Kmer_BankId_Count{
    	Type _type;
    	u_int32_t _bankId;
    	u_int16_t _count;

    	Kmer_BankId_Count(){

    	}

    	Kmer_BankId_Count(Type type, u_int64_t bankId, u_int64_t count){
    		_type = type;
    		_bankId = bankId;
    		_count = count;
    	}
    };

    //typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;
	//struct kxpcomp { bool operator() (kxp l,kxp r) { return (get<0>(r) < get<0>(l)); } } ;

	//string _outputDir;
	string _outputFilename;
	//vector<string>& _datasetToMergeDirs;
	//size_t _partitionId;
	Bag<Kmer_BankId_Count>* _outputGzFile;
	Bag<Kmer_BankId_Count>* _cachedBag;
	map<string, u_int64_t> _dummyIdToOrder;
	//size_t _nbBanks;
	//vector<string> _currentDatasetIds;

	DatasetMergerWriter(size_t partitionId, vector<string>& datasetToMergeDirs, const string& outputDir):
		DiskBasedMergeSort<span>(partitionId, datasetToMergeDirs, _dummyIdToOrder, false)
    {
    	//_outputDir = outputDir;
    	//_partitionId = partitionId;

		_outputFilename = outputDir + "/" + Stringify::format("%i", this->_partitionId) + ".gz";
    	//_outputFilename = _outputDir + "/solid/part_" + Stringify::format("%i", partitionId) + "/__p__" + Stringify::format("%i", mergeId) + ".gz.temp";

    	_outputGzFile = new BagGzFile<Kmer_BankId_Count>(_outputFilename);
    	_cachedBag = new BagCache<Kmer_BankId_Count>(_outputGzFile, 10000);

    	//_nbBanks = 0;
		//_nbBanks = _datasetIds.size();

    	//cout << _outputFilename << endl;
    }

	void process(Type& kmer, u_int64_t bankId, u_int64_t abundance){
		_cachedBag->insert(Kmer_BankId_Count(kmer, bankId, abundance));
	}

	void end(){

		_cachedBag->flush();
    	delete _cachedBag;

    	checkGzFile(_outputFilename);
    	/*
		for(size_t i=0; i<this->_datasetToMergeDirs.size(); i++){
			//cout << _datasetIds[i] << endl;
			//string filename = _datasetToMergeDirs[i];
			string filename = this->_datasetToMergeDirs[i] + "/" + Stringify::format("%i", this->_partitionId) + ".gz";
			//string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			System::file().remove(filename);

			//cout << filename << " " << System::file().doesExist(filename) << endl;
		}

		//cout << _outputFilename << endl;

		string newOutputFilename = this->_datasetToMergeDirs[0] + "/" + Stringify::format("%i", this->_partitionId) + ".gz";
		//newOutputFilename.erase(_outputFilename.size()-5, 5);
		//cout << _outputFilename << "   ->   " << newOutputFilename << endl;
    	System::file().rename(_outputFilename, newOutputFilename); //remove .temp at the end of new merged file
    	_outputFilename = newOutputFilename;
		*/
    	//saveMergeInfos();
    }

	//There is a bug in simka, sometimes a gz file is erroneous at the end
	//It's really rare and I can't find it
	//My bad solution is to read the whole gz file as soon as it is close and a segfault will occur if it has a bad format
	//Of course it's a bad solution because it has a impact on simka performances...
	void checkGzFile(const string& filename){
		IterableGzFile<Kmer_BankId_Count>* gzFile = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);
		Iterator<Kmer_BankId_Count>* it = gzFile->iterator();

		it->first();
		while(!it->isDone()){
			it->next();
		}

		delete it;
		delete gzFile;
	}



};













/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2MergeAlgorithm : public Algorithm
{

public:

    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;

	string _inputFilename;
	string _outputDir;
	size_t _kmerSize;
	size_t _partitionId;
	vector<string> _datasetIds;
	vector<string> _kmerSpectrumDirs;
	size_t _nbCores;
	size_t _nbBanks;
	map<string, string> _links;
	//string _mergeDestDir;
	u_int64_t _nbMergedBanks;
	bool _justSaveMergeInfos;
	string _databaseDir;
	string _tmpDir;

	Simka2MergeAlgorithm(IProperties* options):
		Algorithm("simka", -1, options)
	{
		_nbMergedBanks = 0;
	}

	void execute(){
		parseArgs();
		layoutInputFilename();
		//createIdsOrderIndex();
		//loadLinks();
		merge();
		//partitionKmerCounts();
		//saveLinks();
		//saveInfos();

		//system(("rm -rf " + _outputDirTemp).c_str());
		//System::file().rmdir(_outputDirTemp);

		//cout << "heo" << endl;
		//delete config;
		//cout << "heo" << endl;
		//writeFinishSignal();
	}

	void parseArgs(){

    	_databaseDir =  getInput()->getStr(STR_SIMKA2_DATABASE_DIR);
    	_inputFilename =  getInput()->getStr(STR_URI_INPUT);
    	_partitionId =   getInput()->getInt(STR_SIMKA2_PARTITION_ID);
    	_outputDir =  getInput()->getStr(STR_URI_OUTPUT);
    	_nbCores =  getInput()->getInt(STR_NB_CORES);
    	_justSaveMergeInfos =  getInput()->get("-save-merge-info");
    	//_kmerSize =  getInput()->getStr(STR_URI_INPUT);

    	_tmpDir = _databaseDir + "/merge";
	}


	void layoutInputFilename(){


		//string inputDir = _outputDirTemp + "/input/";
		ifstream inputFile(_inputFilename.c_str());

		//_banksInputFilename =  inputDir + "__input_simka__"; //_inputFilename + "_dsk_dataset_temp__";
		//IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		string line;
		//string linePart;

		//string bankFileContents = "";

		u_int64_t lineIndex = 0;

		while(getline(inputFile, line)){

			line.erase(std::remove(line.begin(),line.end(),' '),line.end());
			if(line == "") continue;

			string datasetID = SimkaIoUtils::getDatasetID(line);//System::file().getBaseName(line);
			//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());


			lineIndex += 1;

			_datasetIds.push_back(datasetID);
			_kmerSpectrumDirs.push_back(line);
			//_bankNames.push_back(bankId);


		}

		_nbBanks = _kmerSpectrumDirs.size();
		//_mergeDestDir = _kmerSpectrumDirs[0];
	}

	/*
	void createIdsOrderIndex(){

		Simka2Database database(_databaseDir);

		for(size_t i=0; i<database._entries.size(); i++){
			_idToOrder[database._entries[i]] = i;
		}
	}*/

	/*
	void loadLinks(){
		for(size_t i=0; i<_nbBanks; i++){

			string datasetID = getDatasetID(_kmerSpectrumDirs[i]);
			//string datasetID = System::file().getBaseName(_kmerSpectrumFilenames[i]);
			//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());

			string mergedLinkinkFilename = _kmerSpectrumDirs[i] + "/merged_link.bin";


			if(System::file().doesExist(mergedLinkinkFilename)){

				//string mergedDir = System::file().getDirectory(mergeFilenames[i]);
				ifstream mergedLinkFile(mergedLinkinkFilename.c_str(), std::ios::binary);
				u_int64_t size;
				mergedLinkFile.read((char*)(&size), sizeof(size));
				std::vector<char> buffer(size); // create a buffer
				mergedLinkFile.read( &buffer[0], buffer.size() ); // read to buffer
				string linkedDatasetID( buffer.begin(), buffer.end() ); // copy from vector
				//mergedLinkFile.read(&datasetID[0], size);
				mergedLinkFile.close();

				_links[datasetID] = linkedDatasetID;
			}
			else{
				_links[datasetID] = datasetID;
			}

			cout << "Link: " << datasetID << " to " << _links[datasetID] << endl;

			//cout << partFilename << " " << getFileSize(partFilename) << endl;
		}
	}*/

	//typedef tuple<u_int64_t, string> sortItem_Size_Filename_ID;

	//static bool sortFileBySize (sortItem_Size_Filename_ID& i, sortItem_Size_Filename_ID& j){
	//	return ( get<0>(i) < get<0>(j) );
	//}



	void merge(){

		if(_justSaveMergeInfos){
			saveMergeInfos();
		}
		else{
			//for(size_t i=0; i<0; i++){
				DatasetMergerWriter<span> diskBasedMergeSort(_partitionId, _kmerSpectrumDirs, _outputDir);
				diskBasedMergeSort.execute();
				//_nbMergedBanks = diskBasedMergeSort._nbBanks;
				//break;
			//}
		}


		//mergeDatasets();
		//for(size_t i=0; i<_nbBanks; i++){
			//cout << i << endl;
			//mergePartition(i);
			//break;
			//exit(1);
		//}

		//_partitionId = p.partitionId;

		//createDatasetIdList(p);
		//_nbBanks = _datasetIds.size();

	}

    void saveMergeInfos(){

    	_nbMergedBanks = 0;

		for(size_t i=0; i<_kmerSpectrumDirs.size(); i++){
			string kmerSpectrumDir = _kmerSpectrumDirs[i];

			u_int64_t nbMergedBanks = 0;
			//size_t bankOffset = 0;
	    	string mergeInfoFilename = kmerSpectrumDir + "/merge_info.bin";

			ifstream mergedLinkFile(mergeInfoFilename.c_str(), std::ios::binary);
			mergedLinkFile.read((char*)(&nbMergedBanks), sizeof(nbMergedBanks));
			mergedLinkFile.close();

			_nbMergedBanks += nbMergedBanks;
		}


		cout << "saving" << endl;
		cout << _nbMergedBanks << endl;






    	//cout << _nbMergedBanks << endl;
    	//cout << "\t\tsaving info:" << endl;
    	//_outputFilename = _datasetToMergeDirs[0] + "/" + Stringify::format("%i", partitionId) + ".temp";


    	string mergeInfoFilenameTemp = _outputDir + "/merge_info.bin";
		ofstream kmerSpectrumInfoFile(mergeInfoFilenameTemp.c_str(), std::ios::binary);

		kmerSpectrumInfoFile.write((char const*)(&_nbMergedBanks), sizeof(_nbMergedBanks));

		for(size_t i=0; i<_kmerSpectrumDirs.size(); i++){
			string kmerSpectrumDir = _kmerSpectrumDirs[i];
			//string filename = kmerSpectrumDir + "/" + Stringify::format("%i", _partitionId) + ".gz";

			//cout << filename << endl;

			u_int64_t nbMergedBanks = 0;
			//size_t bankOffset = 0;
	    	string mergeInfoFilename = kmerSpectrumDir + "/merge_info.bin";

			//cout << kmerSpectrumDir << "  " << System::file().doesExist(mergeInfoFilename)  << "  " << mergeInfoFilename<< endl;


	        //outputInfoFile.write((char const*)(&_nbReads), sizeof(_nbReads));
			//outputInfoFile.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
			//outputInfoFile.write((char const*)(&nbKmers), sizeof(nbKmers));
			//outputInfoFile.write((char const*)(&chord_N2), sizeof(chord_N2));

			//if(System::file().doesExist(mergeInfoFilename)){
			ifstream mergedLinkFile(mergeInfoFilename.c_str(), std::ios::binary);
			mergedLinkFile.read((char*)(&nbMergedBanks), sizeof(nbMergedBanks));
			for(size_t i=0; i<nbMergedBanks; i++){
				//string datasetId = _allDatasetIds[_datasetIds[i]];
				//u_int64_t size = datasetId.size();
				//mergedLinkFile.write((char const*)(&size), sizeof(size));
				//mergedLinkFile.write(datasetId.c_str(), size);
				//cout << datasetId << endl;

				SimkaIoUtils::simka2_transferDatasetInfo(mergedLinkFile, kmerSpectrumInfoFile);

				//cout << size << " " << linkedDatasetID << endl;
				//simka2_writeString(linkedDatasetID, kmerSpectrumInfoFile);
				//_currentDatasetIds.push_back(linkedDatasetID);
				//cout << "\tmerge info:" << linkedDatasetID << endl;
			}
			mergedLinkFile.close();

				//_nbBanks += nbMergedBanks;
				//}
				//else{

	    		//cout << kmerSpectrumDir << endl;
				//string datasetID = getDatasetID(kmerSpectrumDir);

				//cout << datasetID  << endl;
				//simka2_writeString(datasetID, kmerSpectrumInfoFile);

				//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
				//string datasetID = System::file().getBaseName(filename);
				//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
				//_currentDatasetIds.push_back(datasetID);

				//nbMergedBanks = 1;
				//}


		}

		kmerSpectrumInfoFile.close();

    	//string mergeInfoFilename = _outputDir + "/merge_info.bin";
    	//System::file().remove(mergeInfoFilename);
    	//System::file().rename(mergeInfoFilenameTemp, mergeInfoFilename);
		/*
		u_int64_t nbBanks = _kmerSpectrumDirs.size();

		kmerSpectrumInfoFile.write((char const*)(&nbBanks), sizeof(nbBanks));
		for(size_t i=0; i<_kmerSpectrumDirs.size(); i++){
			string datasetId = getDatasetID(_kmerSpectrumDirs[i]);
			u_int64_t size = datasetId.size();
			kmerSpectrumInfoFile.write((char const*)(&size), sizeof(size));
			kmerSpectrumInfoFile.write(datasetId.c_str(), size);
	    	cout << "\t\t" << datasetId << endl;
			//cout << datasetId << endl;
		}
		kmerSpectrumInfoFile.close();
		*/
		//cout << "\tInfo saved:   " << mergeInfoFilename << endl;
    }



	/*
	//void mergePartition(size_t partitionID){
	void mergeDatasets(){

		//vector<string> filenames;
		vector<sortItem_Size_Filename_ID> filenameSizes;
		vector<string>  usedDatasetIds;

		for(size_t i=0; i<_nbBanks; i++){

			string datasetID = getDatasetID(_kmerSpectrumDirs[i]);

			cout << datasetID <<  " " << _links[datasetID] << endl;
			datasetID = _links[datasetID];
			//cout << datasetID << endl;
			if (find (usedDatasetIds.begin(), usedDatasetIds.end(), datasetID) != usedDatasetIds.end()){
				continue;
			}
			usedDatasetIds.push_back(datasetID);

			string linkedFilename = System::file().getDirectory(_kmerSpectrumDirs[i]) + "/" + datasetID;

			//filenames.push_back(filename);
			filenameSizes.push_back(sortItem_Size_Filename_ID(getDatasetSize(linkedFilename), linkedFilename));

			cout << "used: " << datasetID << endl;
		}

		//cout << filenameSizes.size() << endl;
		//exit(1);
		//string partDir = p.outputDir + "/solid/part_" + Stringify::format("%i", _partitionId) + "/";
		//vector<string> filenames = System::file().listdir(partDir);
		//cout << filenames.size() << endl;
		vector<string> partFilenames;

		//cout << "mettre un while ici" << endl;

		size_t nbPartLimit = 2;
		while(filenameSizes.size() > nbPartLimit){

			//cout << "Start merging pass" << endl;
			sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);

			//vector<size_t> mergeDatasetIds;
			vector<string> mergeDatasetDirs;
			//vector<size_t> toRemoveItem;

			cout << endl << "Merging ";
			//cout << endl;
			//cout << "merging" << endl;
			for(size_t i=0; i<nbPartLimit; i++){
				sortItem_Size_Filename_ID sfi = filenameSizes[i];
				//mergeDatasetIds.push_back(get<2>(sfi));
				mergeDatasetDirs.push_back(get<1>(sfi));
				//cout << get<1>(sfi) << endl;
				//datasetIndex += 1;
				//if(datasetIndex >= _nbBanks) break;

				cout << getDatasetID(mergeDatasetDirs[i]) << " ";
				//cout << "First val must never be greater than second:   " << i << "  " << _nbBanks << endl;
				//cout << "\t" << get<1>(sfi) << endl;
			}

			string mergeDestID = getDatasetID(mergeDatasetDirs[0]);
			cout << "    in    " << mergeDestID << endl;


			//size_t mergedId = mergeDatasetIds[0];
			//string mergeOuputDir = mergeDatasetDirs[0];
			//mergeOuputDir = System::file().getDirectory(mergeOuputDir);
			//mergeOuputDir.erase(mergeOuputDir.end()-1);
			//mergeOuputDir += "_merged";

			//System::file().mkdir(mergeOuputDir, -1);

			for(size_t i=0; i<mergeDatasetDirs.size(); i++){

				_links[getDatasetID(mergeDatasetDirs[i])] = mergeDestID;
				updateLinks(getDatasetID(mergeDatasetDirs[i]), mergeDestID);

				filenameSizes.erase(filenameSizes.begin());
			}

			mergePartitions(mergeDatasetDirs);


			//DiskBasedMergeSort<span> diskBasedMergeSort(mergedId, _datasetIds[mergedId], mergeOuputDir, mergeDatasetIds, partitionID);


			//cout << "\t" << mergeOuputDir << endl;

			filenameSizes.push_back(sortItem_Size_Filename_ID(getDatasetSize(mergeDestID), mergeDatasetDirs[0]));

			//exit(1);
		}

	}

	void mergePartitions(vector<string>& mergeDatasetDirs){

		for(size_t i=0; i<SIMKA2_NB_PARTITIONS; i++){
			DiskBasedMergeSort<span> diskBasedMergeSort(i, mergeDatasetDirs);
			diskBasedMergeSort.execute();
			//break;
			exit(1);
		}
	}



    void loadMergeInfos(){

    }

	void updateLinks(const string& oldId, const string&  newId){

		for(size_t i=0; i<_nbBanks; i++){
			if(_links[_datasetIds[i]] == oldId){
				_links[_datasetIds[i]] = newId;
			}
		}
	}

	void saveLinks(){

		//cout << _nbBanks << endl;
		for(size_t i=0; i<_nbBanks; i++){

			cout << _datasetIds[i] << "  " << _links[_datasetIds[i]] << endl;

			string datasetID = _links[_datasetIds[i]];
			//string datasetID = System::file().getBaseName(_kmerSpectrumFilenames[i]);
			//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());

			string mergedLinkinkFilename = _kmerSpectrumDirs[i] + "/merged_link.bin";


			//string mergedDir = System::file().getDirectory(mergeFilenames[i]);
			ofstream mergedLinkFile(mergedLinkinkFilename.c_str(), std::ios::binary);
			//cout << datasetID << endl;
			u_int64_t size = datasetID.size();
			//cout << "-" << size << endl;
			//cout << sizeof(size) << endl;
			mergedLinkFile.write((char const*)(&size), sizeof(size));
			mergedLinkFile.write(datasetID.c_str(), size);
			mergedLinkFile.close();

		}
	}

	void updateInfos(){

	}

	void saveInfos(){

	}*/

	void writeFinishSignal(){
		string finishFilename = _tmpDir + "/" + Stringify::format("%i", _partitionId) + "-success";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
	}



};











class Simka2Merge: public Tool{
public:

	string _execFilename;

	Simka2Merge(string execFilename): Tool ("Simka2-Merge"){
		_execFilename = execFilename;

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for merged kmer spectrum", true));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for merged k-mer spectrums", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DATABASE_DIR, "dir path to a simka database", true));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input filename of k-mer spectrums to merge | TODO SPECIF", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_PARTITION_ID, "number of the partition", true));
	    parser->push_front (new OptionNoParam ("-save-merge-info", ""));

	    IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", true));

	    //parser->getParser(STR_NB_CORES)->setVisible(false);
		parser->push_back(kmerParser);
	    /*
		//Main parser
	    //parser->push_front (new OptionNoParam (STR_SIMKA_COMPUTE_DATA_INFO, "compute (and display) information before running Simka, such as the number of reads per dataset", false));
	    //parser->push_front (new OptionNoParam (STR_SIMKA_KEEP_TMP_FILES, "keep temporary files", false));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for kmer spectrum", false, "./spectrum"));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input filename | TODO SPECIF", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DATASET_ID, "identifier of the input dataset", true));


	    //parser->push_back (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
		//IOptionsParser* parser = getParser();
		//IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();
		//parser->push_back(dskParser);
		//dskParser->setVisible(false);
		//cout << parser->getParser(STR_NB_CORES) << endl;
		//

		//parser->push_back(new OptionOneParam(parser->getParser(STR_NB_CORES)->getName(), parser->getParser(STR_NB_CORES)->getHelp(), false, "0"));
		//parser->push_front(dskParser->getParser (STR_URI_OUTPUT_TMP));
		//dskParser->getParser (STR_URI_OUTPUT_TMP)->setMandatory
	    //parser->push_front(dskParser->getParser (STR_URI_OUTPUT));
	    //parser->getParser (STR_URI_OUTPUT)->setHelp("output directory for result files (similarity matrix, heatmaps)");
	    //parser->push_front(dskParser->getParser (STR_URI_INPUT));
	    //parser->getParser(STR_URI_INPUT)->setHelp("input file of datasets. One dataset per line: id filename1 filename2...");

	    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_URI_OUTPUT_TMP)))  {  p->s; }

		//Distance parser
	    //IOptionsParser* distanceParser = new OptionsParser ("distance");
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES, "compute all simple distances (Chord, Hellinger...)", false));
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES, "compute all complex distances (Jensen-Shannon...)", false));


		//Kmer parser
	    IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));
	    //kmerParser->push_back(dskParser->getParser (STR_KMER_SIZE));
	    //kmerParser->push_back(new OptionOneParam (STR_KMER_PER_READ.c_str(), "number of selected kmers per read", false, "0"));
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "1"));
	    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "2"));
	    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "max abundance a kmer can have to be considered", false, "999999999"));
	    //kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MIN));
	    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
	    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

	    //kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MAX));
	    //kmerParser->push_back(dskParser->getParser (STR_SOLIDITY_KIND));
	    //kmerParser->getParser (STR_SOLIDITY_KIND)->setHelp("TODO");
	    //kmerParser->push_back (new OptionNoParam (STR_SIMKA_SOLIDITY_PER_DATASET.c_str(), "do not take into consideration multi-counting when determining solid kmers", false ));
	    //kmerParser->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX.c_str(), "minimal Shannon index a kmer should have to be kept. Float in [0,2]", false, "0" ));


	    //Read filter parser
	    IOptionsParser* readParser = new OptionsParser ("read");
	    readParser->push_back (new OptionOneParam (STR_SIMKA_MAX_READS.c_str(), "maximum number of reads to process. Set to 0 to use all reads", false, "0" ));
	    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE.c_str(), "minimal size a read should have to be kept", false, "0" ));
	    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(), "minimal Shannon index a read should have to be kept. Float in [0,2]", false, "0" ));

	    //Core parser
	    IOptionsParser* coreParser = new OptionsParser ("core");
	    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
	    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (MB)", false, "8000"));
	    //coreParser->push_back(dskParser->getParser ());
	    //coreParser->push_back(dskParser->getParser (STR_MAX_DISK));

	    //Distances
	    //IOptionsParser* distanceParser = new OptionsParser ("distances");
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_BRAYCURTIS.c_str(), "compute Bray Curtis distance"));
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_CHORD.c_str(), "compute Chord distance"));
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_HELLINGER.c_str(), "compute Hellinger distance"));
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_CANBERRA.c_str(), "compute Canberra distance"));
	    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_KULCZYNSKI.c_str(), "compute Kulczynski distance"));


		//parser->push_back(distanceParser);
		parser->push_back(kmerParser);
		parser->push_back(readParser);
		parser->push_back(coreParser);
		//parser->push_back(distanceParser);


		IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();

	    //if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE)))  {  p->setDefaultValue ("7"); }
		parser->push_back(dskParser);
		dskParser->setVisible(false);

		parser->getParser(STR_NB_CORES)->setVisible(false);
		//getParser()->push_back(parser);
	    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

	    //return parser;*/
	}

	~Simka2Merge(){

	}

	struct Parameter
	{
	    //Parameter (Simka& simka, IProperties* props) : props(props) {}
	    Parameter (IProperties* props, const string& execFilename) : _props(props), _execFilename(execFilename) {}
	    //Simka& _simka;
	    IProperties* _props;
	    string _execFilename;
	};

	template<size_t span> struct Functor {

    	void operator ()  (Parameter p){

    		Simka2MergeAlgorithm<span>* algo = new Simka2MergeAlgorithm<span>(p._props);
    		algo->execute();
    		delete algo;
    	}


	};

	void execute ()
	{
		IProperties* input = getInput();
		//Parameter params(*this, getInput());
		Parameter params(input, _execFilename);

		size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

	    Integer::apply<Functor,Parameter> (kmerSize, params);
	}

};









int main (int argc, char* argv[])
{
    try
    {
    	Simka2Merge(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



