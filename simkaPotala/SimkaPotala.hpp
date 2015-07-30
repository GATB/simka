/*
 * SimkaFusion.hpp
 *
 *  Created on: 22 juil. 2015
 *      Author: gbenoit
 */

#ifndef TOOLS_SIMKA_SRC_SIMKAFUSION_HPP_
#define TOOLS_SIMKA_SRC_SIMKAFUSION_HPP_

#include <gatb/gatb_core.hpp>
#include <SimkaAlgorithm.hpp>
#include <KmerCountCompressor.hpp>
#include <Simka.hpp>

//#define SERIAL
#define SLEEP_TIME_SEC 1


const string STR_SIMKA_JOB_COMMAND_COUNT = "-cmd-job-count";
const string STR_SIMKA_JOB_COMMAND_MERGE = "-cmd-job-merge";
const string STR_SIMKA_NB_JOB_COUNT = "-max-job-count";
const string STR_SIMKA_NB_JOB_MERGE = "-max-job-merge";


class SimkaBankSample : public BankDelegate
{
public:


	SimkaBankSample (IBank* ref) : BankDelegate (ref)  {
	}

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

        Iterator<Sequence>* it = _ref->iterator ();

        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

    	TruncateIterator<Sequence>* truncIt = new TruncateIterator<Sequence>(*iterators[0], 10000);
    	return truncIt;
    }

private:

};







template<size_t span>
class SimkaNullProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

    SimkaNullProcessor(){}
	~SimkaNullProcessor(){}
    CountProcessorAbstract<span>* clone ()  {  return new SimkaNullProcessor ();  }
	void finishClones (vector<ICountProcessor<span>*>& clones){}
	bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum){return false;}


};

template<size_t span>
class SimkaCompProcessor : public CountProcessorAbstract<span>{

public:


    SimkaCompProcessor(KmerCountCompressor<span>* comp){
    	_comp = comp;
    }
	~SimkaCompProcessor(){}
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCompProcessor (_comp);  }
	void finishClones (vector<ICountProcessor<span>*>& clones){}
	bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum){
		_comp->insert(partId, kmer, count);
		return true;
	}

private:
	KmerCountCompressor<span>* _comp;

};



template<size_t span>
class SimkaFusion {
public:


	SimkaFusion(IProperties* options){

		_maxNbReads = 0;
		_options = options;

		_inputFilename = _options->getStr(STR_URI_INPUT);
		_outputDir = _options->getStr(STR_URI_OUTPUT);
		_outputDirTemp = _options->getStr(STR_URI_OUTPUT_TMP);

		_maxJobCount = _options->getInt(STR_SIMKA_NB_JOB_COUNT);
		_maxJobMerge = _options->getInt(STR_SIMKA_NB_JOB_MERGE);
		_cmdJobCount = _options->getStr(STR_SIMKA_JOB_COMMAND_COUNT);
		_cmdJobMerge = _options->getStr(STR_SIMKA_JOB_COMMAND_MERGE);



		//string solidFilename = _outputDir + "/solid/" +  p.bankName + suffix + ".h5";

		//cout << "SimkaFusion constructor       " << _outputDirTempFilter << endl;
	}

	~SimkaFusion(){

	}


	void execute(){


		//System::file().rmdir(_outputDirTemp + "/input/");
		//System::file().rmdir(_outputDirTemp + "/solid/");
		//System::file().rmdir(_outputDir + "/temp/");
		//System::file().rmdir(_outputDir);


		//_stats = new SimkaStatistics(_nbBanks);


		layoutInputFilename_FROM_SIMKA();
		computeMaxReads();
		createDirs();
		layoutInputFilename();


		createConfig();

		sleep(SLEEP_TIME_SEC);

		count();

		sleep(SLEEP_TIME_SEC);

		merge();

		sleep(SLEEP_TIME_SEC);

		stats();

	}

	void layoutInputFilename_FROM_SIMKA(){

		if(_options->getInt(STR_VERBOSE) != 0){
			cout << endl << "Creating input" << endl;
		}

		_banksInputFilename = _inputFilename + "_dsk_dataset_temp__";
		IFile* inputFile = System::file().newFile(_inputFilename, "rb");
		IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		inputFile->seeko(0, SEEK_END);
		u_int64_t size = inputFile->tell();
		inputFile->seeko(0, SEEK_SET);
		char buffer2[size];
		inputFile->fread(buffer2, size, size);
		string fileContents(buffer2, size);

		string line;
		string linePart;
		vector<string> linePartList;
		stringstream fileContentsStream(fileContents);

		string bankFileContents = "";

		u_int64_t lineIndex = 0;

		while(getline(fileContentsStream, line)){

			if(line == "") continue;

			stringstream lineStream(line);
			linePartList.clear();
			//vector<string> filenames;

			while(getline(lineStream, linePart, ' ')){

				if(linePart != ""){
					linePartList.push_back(linePart);
				}
			}

	/*
			bool valid = true;
			//cout << linePartList.size() << endl;
			//Bank id
			for(size_t i=1; i<linePartList.size(); i++){
				string filename = linePartList[1];
				//cout << filename << endl;
				if( ! System::file().doesExist(filename)){
					cout << "\tFilename does not exist: " << filename << endl;
					valid = false;
					//break;
				}
			}

			if(!valid){
				continue;
			}*/

			string bankId = linePartList[0];
			_bankNames.push_back(bankId);


			 //ID and one filename
			if(linePartList.size() == 2){
				bankFileContents += linePartList[1] + "\n";
				_nbBankPerDataset.push_back(1);
			}
			//ID and list of filename (paired files for example)
			else{
				char buffer[200];
				snprintf(buffer,200,"%llu", lineIndex);
				string subBankFilename = _banksInputFilename + "_" + string(buffer);
				_tempFilenamesToDelete.push_back(subBankFilename);
				IFile* subBankFile = System::file().newFile(subBankFilename, "wb");
				string subBankContents = "";

				for(size_t i=1; i<linePartList.size(); i++){
					subBankContents += linePartList[i] + "\n";
				}
				subBankContents.erase(subBankContents.size()-1);
				//subBankContents.pop_back(); // "remove last /n
				subBankFile->fwrite(subBankContents.c_str(), subBankContents.size(), 1);
				subBankFile->flush();
				delete subBankFile;

				bankFileContents += System::file().getBaseName(subBankFilename) + "\n";
				_nbBankPerDataset.push_back(linePartList.size() - 1); //linePartList.size() - 1 = nb sub banks
				//_nbReadsPerDataset.push_back(ceil(_maxNbReads / (float)()));
			}

			lineIndex += 1;
		}

		bankFileContents.erase(bankFileContents.size()-1);
		//bankFileContents.pop_back(); // "remove last /n

		bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);

		delete inputFile;
		bankFile->flush();
		delete bankFile;

		//for(int i=0; i<_nbBanksOfDataset.size(); i++){
		//	cout << i << "   "  << _nbBanksOfDataset[i] << endl;
		//}

		if(_options->getInt(STR_VERBOSE) != 0){
			cout << "\tNb input datasets: " << _bankNames.size() << endl;
		}

		cout << endl;

	}


	void layoutInputFilename(){

		//if(_options->getInt(STR_VERBOSE) != 0){
		//	cout << endl << "Creating input" << endl;
		//}

		IFile* inputFile = System::file().newFile(_inputFilename, "rb");

		//Hold dataset ids
		string datasetIdFilename = _outputDirTempFilter + "/" + "datasetIds";
		IFile* datasetIdFile = System::file().newFile(datasetIdFilename, "wb");

		inputFile->seeko(0, SEEK_END);
		u_int64_t size = inputFile->tell();
		inputFile->seeko(0, SEEK_SET);
		char buffer2[size];
		inputFile->fread(buffer2, size, size);
		string fileContents(buffer2, size);

		string line;
		string linePart;
		vector<string> linePartList;
		stringstream fileContentsStream(fileContents);

		//string bankFileContents = "";

		u_int64_t lineIndex = 0;

		while(getline(fileContentsStream, line)){

			if(line == "") continue;

			stringstream lineStream(line);
			linePartList.clear();
			//vector<string> filenames;

			while(getline(lineStream, linePart, ' ')){

				if(linePart != ""){
					linePartList.push_back(linePart);
				}
			}

	/*
			bool valid = true;
			//cout << linePartList.size() << endl;
			//Bank id
			for(size_t i=1; i<linePartList.size(); i++){
				string filename = linePartList[1];
				//cout << filename << endl;
				if( ! System::file().doesExist(filename)){
					cout << "\tFilename does not exist: " << filename << endl;
					valid = false;
					//break;
				}
			}

			if(!valid){
				continue;
			}*/

			string bankId = linePartList[0];

			string bankIdLine = bankId + '\n';
			datasetIdFile->fwrite(bankIdLine.c_str(), bankIdLine.size(), 1);

			//_bankNames.push_back(bankId);

			IFile* subBankFile = System::file().newFile(_outputDirTempFilter + "/input/" + bankId, "wb");
			string subBankContents = "";

			for(size_t i=1; i<linePartList.size(); i++){
				subBankContents += linePartList[i] + "\n";
			}
			subBankContents.erase(subBankContents.size()-1);
			//subBankContents.pop_back(); // "remove last /n
			subBankFile->fwrite(subBankContents.c_str(), subBankContents.size(), 1);
			subBankFile->flush();
			delete subBankFile;

			 //ID and one filename
			if(linePartList.size() == 2){
				//bankFileContents += linePartList[1] + "\n";
				_nbBankPerDataset.push_back(1);
			}
			//ID and list of filename (paired files for example)
			else{
				//char buffer[200];
				//snprintf(buffer,200,"%llu", lineIndex);
				//string subBankFilename = _banksInputFilename + "_" + string(buffer);
				//_tempFilenamesToDelete.push_back(subBankFilename);


				//bankFileContents += System::file().getBaseName(subBankFilename) + "\n";
				_nbBankPerDataset.push_back(linePartList.size() - 1); //linePartList.size() - 1 = nb sub banks
				//_nbReadsPerDataset.push_back(ceil(_maxNbReads / (float)()));
			}

			lineIndex += 1;
		}

		//bankFileContents.erase(bankFileContents.size()-1);
		//bankFileContents.pop_back(); // "remove last /n

		//bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);
		datasetIdFile->flush();
		delete datasetIdFile;
		delete inputFile;
		//bankFile->flush();
		//delete bankFile;

		//for(int i=0; i<_nbBanksOfDataset.size(); i++){
		//	cout << i << "   "  << _nbBanksOfDataset[i] << endl;
		//}

		//if(_options->getInt(STR_VERBOSE) != 0){
		//	cout << "\tNb input datasets: " << _bankNames.size() << endl;
		//}

		cout << endl;

		//_nbBanks = _bankNames.size();




	}

	void computeMaxReads(){
		IBank* bank = Bank::open(_banksInputFilename);
		LOCAL(bank);
		_nbBanks = bank->getCompositionNb();

		//cout << _banksInputFilename << endl;
		if(_maxNbReads == 0){
			cout << bank->estimateNbItems() << endl;
			if(_options->getInt(STR_VERBOSE) != 0)
				cout << "-maxNbReads is not defined. Simka will estimating it..." << endl;
			_maxNbReads = bank->estimateNbItems() / _nbBanks;
			_maxNbReads -= (_maxNbReads/10);
			if(_options->getInt(STR_VERBOSE) != 0)
				cout << "Max nb reads: " << _maxNbReads << endl << endl;
		}

		for(size_t i=0; i<_nbBankPerDataset.size(); i++){
			//cout << _maxNbReads << " " << _nbBankPerDataset[i] << endl;
			_nbReadsPerDataset.push_back( ceil(_maxNbReads / (float)(_nbBankPerDataset[i])) );
		}

		//_nbReadsPerDataset = nbReadsPerDataset;

		System::file().remove(_banksInputFilename);
	    for(size_t i=0; i<_tempFilenamesToDelete.size(); i++){
	    	System::file().remove(_tempFilenamesToDelete[i]);
	    }
	}

	void createDirs(){

		cout << "la" <<  _maxNbReads << endl;
		string suffix = "";
		suffix += "m" + _options->getStr(STR_SIMKA_MIN_READ_SIZE);
		suffix += "_s" + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX);
		suffix += "_n" + SimkaAlgorithm<>::toString(_maxNbReads);
		_outputDirTempFilter = _outputDirTemp + "/" + suffix + "/";

		System::file().mkdir(_outputDirTemp, -1);
		System::file().mkdir(_outputDirTempFilter, -1);
		System::file().mkdir(_outputDirTempFilter + "/input/", -1);
		System::file().mkdir(_outputDirTempFilter + "/solid/", -1);
		System::file().mkdir(_outputDirTempFilter + "/temp/", -1);
		System::file().mkdir(_outputDirTempFilter + "/count_synchro/", -1);
		System::file().mkdir(_outputDirTempFilter + "/merge_synchro/", -1);
		System::file().mkdir(_outputDirTempFilter + "/stats/", -1);

	}

	void createConfig(){
	    Storage* storage = 0;
        storage = StorageFactory(STORAGE_HDF5).create (_outputDirTempFilter + "/" + "config.h5", true, false);
        LOCAL (storage);

        IBank* bank = Bank::open(_outputDirTempFilter + "/input/" + _bankNames[0]);
        IBank* sampleBank = new SimkaBankSample(bank);
		SortingCountAlgorithm<span> sortingCount (sampleBank, _options);

		SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();

		sortingCount.addProcessor (proc);

		// We launch the algorithm
		sortingCount.execute();

		Configuration config = sortingCount.getConfig();
		config.save(storage->getGroup(""));
		sortingCount.getRepartitor()->save(storage->getGroup(""));

		_nbPartitions = config._nb_partitions;
	    //setStorage (storage);

		//delete storage;
		//sampleBank->forget();
	}


	void count(){

		cout << "Counting datasets..." << endl;

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<_bankNames.size(); i++){

	    	cout << "\tDataset " << i << endl;

			string finishFilename = _outputDirTempFilter + "/count_synchro/" +  _bankNames[i] + ".ok";
			if(System::file().doesExist(finishFilename)){
				cout << "\t" << _bankNames[i] << " already counted (remove file " << finishFilename << " to count again)" << endl;
			}
			else{

				filenameQueue.push_back(_bankNames[i]);


				string tempDir = _outputDirTempFilter + "/temp/" + _bankNames[i];

				string command = "";
#ifndef SERIAL
				//command += "nohup";
#endif
				command += " ./simkaCount ";
				command += " " + string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE);
				command += " " + string("-out-tmp-simka") + " " + _outputDirTempFilter;
				command += " " + string("-out-tmp") + " " + tempDir;
				command += " -bank-name " + _bankNames[i];
				command += " " + string(STR_MAX_MEMORY) + " " + _options->getStr(STR_MAX_MEMORY);
				command += " " + string(STR_NB_CORES) + " " + _options->getStr(STR_NB_CORES);
				command += " " + string(STR_URI_INPUT) + " dummy ";
				command += " " + string(STR_KMER_ABUNDANCE_MIN) + " " + _options->getStr(STR_KMER_ABUNDANCE_MIN);
				command += " " + string(STR_SIMKA_MIN_READ_SIZE) + " " + _options->getStr(STR_SIMKA_MIN_READ_SIZE);
				command += " " + string(STR_SIMKA_MIN_READ_SHANNON_INDEX) + " " + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX);
				command += " " + string(STR_SIMKA_MAX_READS) + " " + SimkaAlgorithm<>::toString(_nbReadsPerDataset[i]);
#ifndef SERIAL
				command += " &";
#endif

				command = _cmdJobCount + " " + command;

				cout << "\t" << command << endl;

				system(command.c_str());
				nbJobs += 1;
			}

			if(nbJobs >= _maxJobCount){
				while(true){
					bool isJobAvailbale = false;

					for(size_t j=0; j<filenameQueue.size(); j++){

						string finishFilename2 = _outputDirTempFilter + "/count_synchro/" + filenameQueue[j] + ".ok";
						if(System::file().doesExist(finishFilename2)){
							filenameQueueToRemove.push_back(filenameQueue[j]);
							isJobAvailbale = true;
							nbJobs -= 1;
						}
					}

					if(isJobAvailbale){
						for(size_t j=0; j<filenameQueueToRemove.size(); j++){
							filenameQueue.erase(std::remove(filenameQueue.begin(), filenameQueue.end(), filenameQueueToRemove[j]), filenameQueue.end());
						}
						filenameQueueToRemove.clear();
						break;
					}
					else{
						sleep(1);
					}

					if(i >= _bankNames.size()) break;
				}
			}
	    }

	    cout << nbJobs << endl;

	    while(nbJobs > 0){
			bool isJobAvailbale = false;

			for(size_t j=0; j<filenameQueue.size(); j++){

				string finishFilename2 = _outputDirTempFilter + "/count_synchro/" + filenameQueue[j] + ".ok";
				if(System::file().doesExist(finishFilename2)){
					filenameQueueToRemove.push_back(filenameQueue[j]);
					isJobAvailbale = true;
					nbJobs -= 1;
				}
			}

			if(isJobAvailbale){
				for(size_t j=0; j<filenameQueueToRemove.size(); j++){
					filenameQueue.erase(std::remove(filenameQueue.begin(), filenameQueue.end(), filenameQueueToRemove[j]), filenameQueue.end());
				}
				filenameQueueToRemove.clear();
				//break;
			}
			else{
				sleep(1);
			}
	    }


	    cout << nbJobs << endl;

	}

	void waitForJob(){

	}

	void merge(){
		cout << "Merging datasets..." << endl;

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<_nbPartitions; i++){

	    	cout << "\tDataset " << i << endl;

	    	string datasetId = SimkaAlgorithm<>::toString(i);
			string finishFilename = _outputDirTempFilter + "/merge_synchro/" +  datasetId + ".ok";
			if(System::file().doesExist(finishFilename)){
				System::file().remove(finishFilename);
			//	cout << "\t" << _bankNames[i] << " already  (remove file " << finishFilename << " to count again)" << endl;
			}
			//else{

				filenameQueue.push_back(datasetId);

				string command = "";

#ifndef SERIAL
				//command += "nohup";
#endif
				command += " ./simkaMerge ";
				command += " " + string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE);
				command += " " + string(STR_URI_INPUT) + " " + _inputFilename;
				command += " " + string("-out-tmp-simka") + " " + _outputDirTempFilter;
				command += " -partition-id " + SimkaAlgorithm<>::toString(i);
				command += " " + string(STR_MAX_MEMORY) + " " + _options->getStr(STR_MAX_MEMORY);
				command += " " + string(STR_NB_CORES) + " " + _options->getStr(STR_NB_CORES);
#ifndef SERIAL
				command += " &";
#endif

				command = _cmdJobMerge + " " + command;

				cout << "\t" << command << endl;

				system(command.c_str());
				nbJobs += 1;
			//}

			if(nbJobs >= _maxJobMerge){


				while(true){

					//cout << nbJobs  << " " << filenameQueue.size()<< endl;
					//cout << "allo..???" << endl;
					bool isJobAvailbale = false;

					for(size_t j=0; j<filenameQueue.size(); j++){

						string finishFilename2 = _outputDirTempFilter + "/merge_synchro/" + filenameQueue[j] + ".ok";
						if(System::file().doesExist(finishFilename2)){
							filenameQueueToRemove.push_back(filenameQueue[j]);
							isJobAvailbale = true;
						}
					}

					if(isJobAvailbale){
						for(size_t j=0; j<filenameQueueToRemove.size(); j++){
							nbJobs -= 1;
							filenameQueue.erase(std::remove(filenameQueue.begin(), filenameQueue.end(), filenameQueueToRemove[j]), filenameQueue.end());
						}
						filenameQueueToRemove.clear();
						break;
					}
					else{
						sleep(1);
					}

					if(i >= _nbPartitions){break;}
				}
			}
	    }

	    cout << nbJobs << endl;

	    while(nbJobs > 0){
			bool isJobAvailbale = false;

			for(size_t j=0; j<filenameQueue.size(); j++){

				string finishFilename2 = _outputDirTempFilter + "/merge_synchro/" + filenameQueue[j] + ".ok";
				if(System::file().doesExist(finishFilename2)){
					filenameQueueToRemove.push_back(filenameQueue[j]);
					isJobAvailbale = true;
				}
			}

			if(isJobAvailbale){
				for(size_t j=0; j<filenameQueueToRemove.size(); j++){
					nbJobs -= 1;
					filenameQueue.erase(std::remove(filenameQueue.begin(), filenameQueue.end(), filenameQueueToRemove[j]), filenameQueue.end());
				}
				filenameQueueToRemove.clear();
			}
			else{
				sleep(1);
			}
	    }


	    cout << nbJobs << endl;

	}

	void stats(){
		cout << endl << "Computing stats..." << endl;
		cout << _nbBanks << endl;

		SimkaStatistics mainStats(_nbBanks);

		for(size_t i=0; i<_nbPartitions; i++){

			Storage* storage = StorageFactory(STORAGE_HDF5).load (_outputDirTempFilter + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".stats");
			LOCAL (storage);

			SimkaStatistics stats(_nbBanks);
			stats.load(storage->getGroup(""));

			cout << stats._nbKmers << endl;
			mainStats += stats;
		}

		mainStats.print();

		mainStats.outputMatrix(_outputDir, _bankNames);
	}

	IBank* _banks;
	IProperties* _options;
	string _inputFilename;
	string _outputDir;
	string _outputDirTemp;
	string _outputDirTempFilter;
	vector<string> _bankNames;
    vector<size_t> _nbBankPerDataset;
    size_t _nbPartitions;
    size_t _nbBanks;
	vector<u_int64_t> _nbReadsPerDataset;
	string _banksInputFilename;
	vector<string> _tempFilenamesToDelete;
	u_int64_t _maxNbReads;
	//IBank* _sampleBank;

	size_t _maxJobCount;
	size_t _maxJobMerge;
	string _cmdJobCount;
	string _cmdJobMerge;

};














/*
IOptionsParser* Simka::createOptionsParser (IOptionsParser* parent)
{
    IOptionsParser* parser = parent; //new OptionsParser ("Simka");

	//IOptionsParser* parser = getParser();
	IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();
	parser->push_back(dskParser);
	dskParser->setVisible(false);
	//cout << parser->getParser(STR_NB_CORES) << endl;
	parser->getParser(STR_NB_CORES)->setVisible(false);

	//Main parser
	parser->push_front(dskParser->getParser (STR_URI_OUTPUT_TMP));
	parser->push_front(dskParser->getParser (STR_URI_OUTPUT));
	parser->getParser (STR_URI_OUTPUT)->setHelp("output directory for result files (similarity matrix, heatmaps)");
	parser->push_front(dskParser->getParser (STR_URI_INPUT));
	parser->getParser(STR_URI_INPUT)->setHelp("input file of datasets. One dataset per line: id filename1 filename2...");


	//Kmer parser
    IOptionsParser* kmerParser = new OptionsParser ("kmer");
    kmerParser->push_back(dskParser->getParser (STR_KMER_SIZE));
    kmerParser->push_back(new OptionOneParam (STR_KMER_PER_READ.c_str(), "number of selected kmers per read", false, "0"));
    kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MIN));
    if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MAX));
    kmerParser->push_back(dskParser->getParser (STR_SOLIDITY_KIND));
    kmerParser->getParser (STR_SOLIDITY_KIND)->setHelp("TODO");
    kmerParser->push_back (new OptionNoParam (STR_SIMKA_SOLIDITY_PER_DATASET.c_str(), "do not take into consideration multi-counting when determining solid kmers", false ));
    kmerParser->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX.c_str(), "minimal Shannon index a kmer should have to be kept. Float in [0,2]", false, "0" ));


    //Read filter parser
    IOptionsParser* readParser = new OptionsParser ("read");
    readParser->push_back (new OptionOneParam (STR_SIMKA_MAX_READS.c_str(), "maximum number of reads per dataset to process", false, "0" ));
    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE.c_str(), "minimal size a read should have to be kept", false, "0" ));
    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(), "minimal Shannon index a read should have to be kept. Float in [0,2]", false, "0" ));

    //Core parser
    IOptionsParser* coreParser = new OptionsParser ("core");
    coreParser->push_back(new OptionOneParam(parser->getParser(STR_NB_CORES)->getName(), parser->getParser(STR_NB_CORES)->getHelp(), false, "0"));
    coreParser->push_back(dskParser->getParser (STR_MAX_MEMORY));
    coreParser->push_back(dskParser->getParser (STR_MAX_DISK));

	parser->push_back(kmerParser);
	parser->push_back(readParser);
	parser->push_back(coreParser);


    return parser;
}*/


class SimkaPotala : public Tool{

public:

	SimkaPotala();
	void execute();
};





#endif
