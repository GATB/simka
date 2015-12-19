/*****************************************************************************
 *   Simka: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2015  INRIA
 *   Authors: G.Benoit, C.Lemaitre, P.Peterlongo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/


#ifndef TOOLS_SIMKA_SRC_SIMKAFUSION_HPP_
#define TOOLS_SIMKA_SRC_SIMKAFUSION_HPP_

#include <gatb/gatb_core.hpp>
#include <SimkaAlgorithm.hpp>
#include <KmerCountCompressor.hpp>
#include <Simka.hpp>

#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>

#include <unistd.h>
#include <sys/wait.h>

//#define CLUSTER
//#define SERIAL
#define SLEEP_TIME_SEC 1


const string STR_SIMKA_CLUSTER_MODE = "-cluster";
const string STR_SIMKA_NB_JOB_COUNT = "-max-count";
const string STR_SIMKA_NB_JOB_MERGE = "-max-merge";
const string STR_SIMKA_JOB_COUNT_COMMAND = "-count-cmd";
const string STR_SIMKA_JOB_MERGE_COMMAND = "-merge-cmd";
const string STR_SIMKA_JOB_COUNT_FILENAME = "-count-file";
const string STR_SIMKA_JOB_MERGE_FILENAME = "-merge-file";

class SimkaBankSample : public BankDelegate
{
public:


	SimkaBankSample (IBank* ref, u_int64_t nbRead) : BankDelegate (ref)  {
		_nbRead = nbRead;
	}

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

        Iterator<Sequence>* it = _ref->iterator ();

        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

    	TruncateIterator<Sequence>* truncIt = new TruncateIterator<Sequence>(*iterators[0], _nbRead);
    	return truncIt;
    }

private:
    u_int64_t _nbRead;
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
class SimkaPotaraAlgorithm : public SimkaAlgorithm<span>{
public:


    typedef typename Kmer<span>::Type  Type;

	SimkaPotaraAlgorithm(IProperties* options):
		SimkaAlgorithm<span>(options)
	{

		_isClusterMode = false;
		//_options = options;

		//_inputFilename = _options->getStr(STR_URI_INPUT);
		//_outputDir = _options->getStr(STR_URI_OUTPUT);
		//_outputDirTemp = _options->getStr(STR_URI_OUTPUT_TMP);

		//_maxNbReads = _options->getInt(STR_SIMKA_MAX_READS);
		//_maxJobCount = _options->getInt(STR_SIMKA_NB_JOB_COUNT);
		//_maxJobMerge = _options->getInt(STR_SIMKA_NB_JOB_MERGE);
		//_jobCountFilename = _options->getStr(STR_SIMKA_JOB_COUNT_FILENAME);
		//_jobMergeFilename = _options->getStr(STR_SIMKA_JOB_MERGE_FILENAME);
		//_jobCountCommand = _options->getStr(STR_SIMKA_JOB_COUNT_COMMAND);
		//_jobMergeCommand = _options->getStr(STR_SIMKA_JOB_MERGE_COMMAND);


		//string solidFilename = _outputDir + "/solid/" +  p.bankName + suffix + ".h5";

		//cout << "SimkaFusion constructor       " << _outputDirTemp << endl;




	}

	~SimkaPotaraAlgorithm(){

	}


	void execute(){

		setup();

		//System::file().rmdir(_outputDirTemp + "/input/");
		//System::file().rmdir(_outputDirTemp + "/solid/");
		//System::file().rmdir(_outputDir + "/temp/");
		//System::file().rmdir(_outputDir);


		//_stats = new SimkaStatistics(_nbBanks);


		//createDirs();
		//layoutInputFilename(_outputDirTemp + "/input/");
		//layoutInputFilename();

		//return;
		//if(!computeMaxReads()) return;

		//sleep(SLEEP_TIME_SEC);

		count();

		//sleep(SLEEP_TIME_SEC);

		merge();

		//sleep(SLEEP_TIME_SEC);

		stats();

	}

	void setup(){
		SimkaAlgorithm<span>::setup();
		parseArgs();
		createDirs();
		layoutInputFilename();
		createConfig();
	}

	void parseArgs() {


		if(this->_options->get(STR_SIMKA_CLUSTER_MODE)){

			cout << "cluster mode activated" << endl;
			cout << "\t-max-memory = memory per job" << endl;
			cout << "\t-nb-cores = cores per job" << endl;
			cout << endl;

			_isClusterMode = true;
			_maxJobCount = this->_options->getInt(STR_SIMKA_NB_JOB_COUNT);
			_maxJobMerge = this->_options->getInt(STR_SIMKA_NB_JOB_MERGE);
			_jobCountFilename = this->_options->getStr(STR_SIMKA_JOB_COUNT_FILENAME);
			_jobMergeFilename = this->_options->getStr(STR_SIMKA_JOB_MERGE_FILENAME);
			_jobCountCommand = this->_options->getStr(STR_SIMKA_JOB_COUNT_COMMAND);
			_jobMergeCommand = this->_options->getStr(STR_SIMKA_JOB_MERGE_COMMAND);

			_maxJobMerge = max((int)_maxJobMerge, (int)30);

			_coresPerJob = this->_nbCores;
			_memoryPerJob = this->_maxMemory;


			IFile* inputFile = System::file().newFile(_jobCountFilename, "rb");
			inputFile->seeko(0, SEEK_END);
			u_int64_t size = inputFile->tell();
			inputFile->seeko(0, SEEK_SET);
			char buffer2[size];
			inputFile->fread(buffer2, size, size);
			string fileContents(buffer2, size);
			_jobCountContents = fileContents;
			delete inputFile;

			inputFile = System::file().newFile(_jobMergeFilename, "rb");
			inputFile->seeko(0, SEEK_END);
			size = inputFile->tell();
			inputFile->seeko(0, SEEK_SET);
			char buffer3[size];
			inputFile->fread(buffer3, size, size);
			string fileContents2(buffer3, size);
			_jobMergeContents = fileContents2;
			delete inputFile;


			return;
		}


		size_t maxCores = System::info().getNbCores();

		_maxJobMerge = maxCores-1;
		size_t maxCoreCount = maxCores-1;

		size_t nbCoresCount = min(maxCoreCount, this->_nbCores);

		u_int64_t minMemory = 2000;
		size_t maxJobCountTemp = this->_maxMemory/minMemory;
		_maxJobCount = min(nbCoresCount, maxJobCountTemp);
		_memoryPerJob = this->_maxMemory / _maxJobCount;
		_coresPerJob = ceil(nbCoresCount / (float)_maxJobCount);


		cout << "Nb jobs count in parallel: " << _maxJobCount << endl;
		cout << "\tCores per jobs: " << _coresPerJob << endl;
		cout << "\tMemory per jobs: " << _memoryPerJob << endl;
		cout << "Nb jobs merge in parallel: " << _maxJobMerge << endl;
		cout << endl;
	}



	void layoutInputFilename(){

		//SimkaAlgorithm<span>::layoutInputFilename();

		string datasetIdFilename = this->_outputDirTemp + "/" + "datasetIds";
		IFile* datasetIdFile = System::file().newFile(datasetIdFilename, "wb");

		for(size_t i=0; i<this->_bankNames.size(); i++){
			string bankName = this->_bankNames[i];

			string bankIdLine = bankName + '\n';
			datasetIdFile->fwrite(bankIdLine.c_str(), bankIdLine.size(), 1);
		}

		datasetIdFile->flush();
		delete datasetIdFile;

	}




	void createDirs(){

		/*
		string suffix = "";
		suffix += "m" + _options->getStr(STR_SIMKA_MIN_READ_SIZE);
		suffix += "_s" + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX);
		suffix += "_n" + SimkaAlgorithm<>::toString(_maxNbReads);
		suffix += "_p" + SimkaAlgorithm<>::toString(_nbAskedPartitions);*/
		//_outputDirTemp = _outputDirTemp; // + "/" + suffix + "/";



		//System::file().mkdir(_outputDirTemp, -1);
		System::file().mkdir(this->_outputDirTemp + "/solid/", -1);
		System::file().mkdir(this->_outputDirTemp + "/temp/", -1);
		System::file().mkdir(this->_outputDirTemp + "/count_synchro/", -1);
		System::file().mkdir(this->_outputDirTemp + "/merge_synchro/", -1);
		System::file().mkdir(this->_outputDirTemp + "/stats/", -1);
		System::file().mkdir(this->_outputDirTemp + "/job_count/", -1);
		System::file().mkdir(this->_outputDirTemp + "/job_merge/", -1);

	}

	void createConfig(){

		string filename = this->_outputDirTemp + "/" + "config.h5";
		if(System::file().doesExist(filename)){

		    try{
				cout << "\tconfig already exists (remove file " << filename << " to config again)" << endl;

				Storage* storage = StorageFactory(STORAGE_HDF5).load (filename);
				LOCAL (storage);
				Configuration* config = new Configuration();
				config->load(storage->getGroup(""));
				_nbPartitions = config->_nb_partitions;
				delete config;

				Repartitor* repartitor = new Repartitor();
				//LOCAL(repartitor);
				repartitor->load(storage->getGroup(""));
				delete repartitor;

				return;
		    }
		    catch (Exception& e)
		    {
		    	cout << "\tcan't open config, computing it again" << endl;
		    	System::file().remove(filename);
		    	createConfig();
		        return;
		    }
		}

		this->_options->setInt(STR_MAX_MEMORY, _memoryPerJob);
		this->_options->setInt(STR_NB_CORES, _coresPerJob);

	    Storage* storage = 0;
        storage = StorageFactory(STORAGE_HDF5).create (filename, true, false);
        LOCAL (storage);


		IBank* bank = Bank::open(this->_outputDirTemp + "/input/" + this->_smallerBankId);
		LOCAL(bank);
        //IBank* bank = Bank::open(_outputDirTemp + "/input/" + _bankNames[0]);
        bank->finalize();
        IBank* sampleBank = new SimkaBankSample(bank, bank->estimateNbItems()/3);
		SortingCountAlgorithm<span> sortingCount (sampleBank, this->_options);

		SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();

		sortingCount.addProcessor (proc);

		// We launch the algorithm
		sortingCount.execute();

		Configuration config = sortingCount.getConfig();


		if(_isClusterMode){
			//config._nb_cached_items_per_core_per_part = 100000;
			_nbPartitions = _maxJobMerge;
			config._nb_partitions = _nbPartitions;

		    uint64_t memoryUsageCachedItems;
		    config._nb_cached_items_per_core_per_part = 1 << 8; // cache at least 256 items (128 here, then * 2 in the next while loop)
		    do
		    {
		        config._nb_cached_items_per_core_per_part *= 2;
		        memoryUsageCachedItems = 1LL * config._nb_cached_items_per_core_per_part *config._nb_partitions * config._nbCores * sizeof(Type);
		    }
		    while (memoryUsageCachedItems < config._max_memory * MBYTE / 10);

		    //cout << config._nb_cached_items_per_core_per_part << endl;
		}
		else{
			_nbPartitions = config._nb_partitions;
		}


        RepartitorAlgorithm<span> repart (bank, storage->getGroup(""), config);
        repart.execute ();
        //setRepartitor (new Repartitor(storage->getGroup("minimizers")));
		//SortingCountAlgorithm<span> sortingCount (sampleBank, _options);




		config.save(storage->getGroup(""));
		//sortingCount.getRepartitor()->save(storage->getGroup(""));
		//delete sampleBank;

	    //setStorage (storage);

		//delete storage;
		//sampleBank->forget();
	}


	void count(){

		vector<string> commands;

		_progress = new ProgressSynchro (
			this->createIteratorListener (this->_bankNames.size(), "Counting datasets"),
			System::thread().newSynchronizer());
		_progress->init ();

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<this->_bankNames.size(); i++){

			string finishFilename = this->_outputDirTemp + "/count_synchro/" +  this->_bankNames[i] + ".ok";
			if(System::file().doesExist(finishFilename)){
				cout << "\t" << this->_bankNames[i] << " already counted (remove file " << finishFilename << " to count again)" << endl;
				continue;
			}
			//else{

			string tempDir = this->_outputDirTemp + "/temp/" + this->_bankNames[i];

			string command = "./simkaCount ";
			command += " " + string(STR_KMER_SIZE) + " " + SimkaAlgorithm<>::toString(this->_kmerSize);
			command += " " + string("-out-tmp-simka") + " " + this->_outputDirTemp;
			command += " " + string("-out-tmp") + " " + tempDir;
			command += " -bank-name " + this->_bankNames[i];
			command += " -nb-datasets " + SimkaAlgorithm<>::toString(this->_nbBankPerDataset[i]);
			command += " " + string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(_memoryPerJob);
			command += " " + string(STR_NB_CORES) + " " + SimkaAlgorithm<>::toString(_coresPerJob);
			command += " " + string(STR_URI_INPUT) + " dummy ";
			command += " " + string(STR_KMER_ABUNDANCE_MIN) + " " + SimkaAlgorithm<>::toString(this->_abundanceThreshold.first);
			command += " " + string(STR_SIMKA_MIN_READ_SIZE) + " " + SimkaAlgorithm<>::toString(this->_minReadSize);
			command += " " + string(STR_SIMKA_MIN_READ_SHANNON_INDEX) + " " + Stringify::format("%f", this->_minReadShannonIndex);
			command += " " + string(STR_SIMKA_MAX_READS) + " " + SimkaAlgorithm<>::toString(this->_maxNbReads);
			//command += " -verbose 0";


			filenameQueue.push_back(this->_bankNames[i]);



			if(_isClusterMode){
				string jobFilename = this->_outputDirTemp + "/job_count/job_count_" + SimkaAlgorithm<>::toString(i) + ".bash";
				IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
				string jobCommand = _jobCountContents + '\n' + '\n';
				jobCommand += command;

				cout << "\t" << jobCommand << endl;

				jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
				jobFile->flush();
				string submitCommand = _jobCountCommand + " " + jobFile->getPath();
				delete jobFile;
				system(submitCommand.c_str());
			}
			else{
				command += " &";
				cout << command << endl;
				system(command.c_str());
			}

			nbJobs += 1;
			//cout << "job started" << endl;

			if(nbJobs >= _maxJobCount){
				while(true){
					bool isJobAvailbale = false;

					for(size_t j=0; j<filenameQueue.size(); j++){

						string finishFilename2 = this->_outputDirTemp + "/count_synchro/" + filenameQueue[j] + ".ok";
						if(System::file().doesExist(finishFilename2)){
							filenameQueueToRemove.push_back(filenameQueue[j]);
							isJobAvailbale = true;
							nbJobs -= 1;
							//cout << "job finished" << endl;
							_progress->inc(1);
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

					if(i >= this->_bankNames.size()) break;
				}
			}


	    }

	    while(nbJobs > 0){
			bool isJobAvailbale = false;

			for(size_t j=0; j<filenameQueue.size(); j++){

				string finishFilename2 = this->_outputDirTemp + "/count_synchro/" + filenameQueue[j] + ".ok";
				if(System::file().doesExist(finishFilename2)){
					filenameQueueToRemove.push_back(filenameQueue[j]);
					isJobAvailbale = true;
					nbJobs -= 1;
					_progress->inc(1);
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

	    _progress->finish();
	    delete _progress;
	}


	void merge(){

		_progress = new ProgressSynchro (
			this->createIteratorListener (_nbPartitions, "Merging datasets"),
			System::thread().newSynchronizer());
		_progress->init ();

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<_nbPartitions; i++){

	    	string datasetId = SimkaAlgorithm<>::toString(i);
			string finishFilename = this->_outputDirTemp + "/merge_synchro/" +  datasetId + ".ok";

			if(System::file().doesExist(finishFilename)){
				cout << "\t" << datasetId << " already merged (remove file " << finishFilename << " to merge again)" << endl;
			}
			else{
				//if(System::file().doesExist(finishFilename)){
				//	System::file().remove(finishFilename);
				//	cout << "\t" << _bankNames[i] << " already  (remove file " << finishFilename << " to count again)" << endl;
				//}
				//else{

				filenameQueue.push_back(datasetId);

				string command = "./simkaMerge ";
				command += " " + string(STR_KMER_SIZE) + " " + SimkaAlgorithm<>::toString(this->_kmerSize);
				command += " " + string(STR_URI_INPUT) + " " + this->_inputFilename;
				command += " " + string("-out-tmp-simka") + " " + this->_outputDirTemp;
				command += " -partition-id " + SimkaAlgorithm<>::toString(i);
				command += " " + string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(this->_maxMemory / this->_nbCores);
				command += " " + string(STR_NB_CORES) + " 1";
				command += " " + string(STR_SIMKA_MIN_KMER_SHANNON_INDEX) + " " + Stringify::format("%f", this->_minKmerShannonIndex);

				SimkaDistanceParam distanceParams(this->_options);
				if(distanceParams._computeBrayCurtis) command += " " + STR_SIMKA_DISTANCE_BRAYCURTIS + " ";
				if(distanceParams._computeCanberra) command += " " + STR_SIMKA_DISTANCE_CANBERRA + " ";
				if(distanceParams._computeChord) command += " " + STR_SIMKA_DISTANCE_CHORD + " ";
				if(distanceParams._computeHellinger) command += " " + STR_SIMKA_DISTANCE_HELLINGER + " ";
				if(distanceParams._computeKulczynski) command += " " + STR_SIMKA_DISTANCE_KULCZYNSKI + " ";

				if(_isClusterMode){
					string jobFilename = this->_outputDirTemp + "/job_merge/job_merge_" + SimkaAlgorithm<>::toString(i) + ".bash";
					IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
					string jobCommand = _jobMergeContents + '\n' + '\n';
					jobCommand += command;

					cout << "\t" << jobCommand << endl;

					jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
					jobFile->flush();
					string submitCommand = _jobMergeCommand + " " + jobFile->getPath();
					delete jobFile;
					system(submitCommand.c_str());
				}
				else{
					command += " &";
					cout << "\t" << command << endl;
					system(command.c_str());
				}

				nbJobs += 1;
			}

			if(nbJobs >= _maxJobMerge){
				while(true){
					bool isJobAvailbale = false;

					for(size_t j=0; j<filenameQueue.size(); j++){

						string finishFilename2 = this->_outputDirTemp + "/merge_synchro/" + filenameQueue[j] + ".ok";
						if(System::file().doesExist(finishFilename2)){
							filenameQueueToRemove.push_back(filenameQueue[j]);
							isJobAvailbale = true;
							nbJobs -= 1;
							_progress->inc(1);
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

					if(i >= this->_bankNames.size()) break;
				}
			}
	    }

	    cout << nbJobs << endl;

	    while(nbJobs > 0){
			bool isJobAvailbale = false;

			for(size_t j=0; j<filenameQueue.size(); j++){

				string finishFilename2 = this->_outputDirTemp + "/merge_synchro/" + filenameQueue[j] + ".ok";
				if(System::file().doesExist(finishFilename2)){
					filenameQueueToRemove.push_back(filenameQueue[j]);
					isJobAvailbale = true;
					nbJobs -= 1;
					_progress->inc(1);
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

	    _progress->finish();
	    delete _progress;
	}

	void stats(){
		cout << endl << "Computing stats..." << endl;
		cout << this->_nbBanks << endl;


		SimkaDistanceParam distanceParams(this->_options);
		SimkaStatistics mainStats(this->_nbBanks, distanceParams);

		for(size_t i=0; i<_nbPartitions; i++){

			Storage* storage = StorageFactory(STORAGE_HDF5).load (this->_outputDirTemp + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".stats");
			LOCAL (storage);

			SimkaStatistics stats(this->_nbBanks, distanceParams);
			stats.load(storage->getGroup(""));

			cout << stats._nbDistinctKmers << endl;
			mainStats += stats;
		}

		mainStats.print();

		mainStats.outputMatrix(this->_outputDir, this->_bankNames);
	}


	//u_int64_t _maxMemory;
	//size_t _nbCores;
	u_int64_t _memoryPerJob;
	size_t _coresPerJob;

	//IBank* _banks;
	//IProperties* _options;
	//string _inputFilename;
	//vector<string> _bankNames;
	//vector<size_t> _nbBankPerDataset;
    size_t _nbPartitions;
    //size_t _nbBanks;
	//vector<u_int64_t> _nbReadsPerDataset;
    //string _banksInputFilename;
    //vector<string> _tempFilenamesToDelete;
    //u_int64_t _maxNbReads;
	//IBank* _sampleBank;

    bool _isClusterMode;
	size_t _maxJobCount;
	size_t _maxJobMerge;
	string _jobCountFilename;
	string _jobMergeFilename;
	string _jobCountCommand;
	string _jobMergeCommand;
	//u_int64_t _nbAskedPartitions;

	string _jobCountContents;
	string _jobMergeContents;

	IteratorListener* _progress;
};












class SimkaPotara : public Tool{

public:

	SimkaPotara();
	void execute();
};

#endif
