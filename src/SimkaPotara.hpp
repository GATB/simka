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
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>

//#include <unistd.h>
//#include <sys/wait.h>
#include <cstdlib>

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



class SimkaBankTemp : public BankDelegate
{
public:

	u_int64_t _refNbReads;
	u_int64_t _refTotalSeqSize;
	u_int64_t _refMaxReadSize;

	/** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
	SimkaBankTemp (IBank* ref, u_int64_t maxReads) : BankDelegate (ref) {

		_maxReads = maxReads;
		//_nbBanks = ref->getCompositionNb();
		ref->estimate(_refNbReads, _refTotalSeqSize, _refMaxReadSize);


    	//cout << _refNbReads << endl;
		//cout << _refTotalSeqSize << endl;
		//cout << _refMaxReadSize << endl;
	}


    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize){


    	if(_maxReads == 0){
    		number = _refNbReads;
    		totalSize = _refTotalSeqSize;
    		maxSize = _refMaxReadSize;
    	}
    	else{

    		u_int64_t maxReads = _maxReads;
    		//u_int64_t maxReads = 0;
    		//for(size_t i=0; i<_nbBanks; i++){
    		//	maxReads += _maxReads * _nbPaireds[i];
    		//}
    		//cout << _refNbReads << endl;
    		//cout << _maxReads*_nbBanks << endl;
    		maxReads = min (maxReads, _refNbReads);
			//cout << "ha " <<  maxReads << endl;

			if(maxReads == _refNbReads){
	    		number = _refNbReads;
	    		totalSize = _refTotalSeqSize;
	    		maxSize = _refMaxReadSize;
			}
			else{
				number = maxReads;
				double factor =  (double)maxReads / (double)_refNbReads;
				totalSize = _refTotalSeqSize * factor;
				maxSize = _refMaxReadSize;
			}
    	}

    	//number = _maxReads;
    	//totalSize = (_totalSizeRef*_nbReadToProcess)/_numberRef;
    	//maxSize = _maxSizeRef;

    	//cout << number2 << endl;

    	//u_int64_t readSize = totalSize2 / number2;
    	//cout << "lal:" << number2 << endl;
    	//number = _maxReads;

    	//number = _nbReadToProcess;
    	//totalSize = _nbReadToProcess*readSize;
    	//maxSize = readSize;

    	//cout << number << endl;
    	//cout << totalSize << endl;
    	//cout << maxSize << endl;
    }

    u_int64_t _maxReads;
};


template<size_t span>
class SimkaPotaraAlgorithm : public SimkaAlgorithm<span>{
public:


    typedef typename Kmer<span>::Type  Type;

	SimkaPotaraAlgorithm(IProperties* options, const string& execFilename):
		SimkaAlgorithm<span>(options)
	{

		_isClusterMode = false;

		//cout << "lala" << endl;
		//cout << _execDir << endl;
       
		if (execFilename == "simka") // i.e. simka binary is found in $PATH
		{
			//cout << endl << "Debug info : execFilename = simka (i.e. simka binary is found in $PATH)" << endl;
			_execDir = "";
		}
		else
		{
			//cout << endl << "Debug info : execFilename = " << execFilename << endl;
			_execDir = System::file().getRealPath(execFilename);
			_execDir = System::file().getDirectory(_execDir) + "/";
		}

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

		parseArgs();

		setup();

		if(!SimkaAlgorithm<span>::isInputValid()) exit(1);

		SimkaAlgorithm<span>::computeMaxReads();

		createConfig();


		count();

		printCountInfo();

		merge();

		stats();


		if(this->_options->getInt(STR_VERBOSE) != 0){
			cout << endl;
			cout << "Output dir: " << this->_outputDir << endl;
			cout << endl;
		}

		//bool keepTempFiles = false;
		if(!this->_keepTmpFiles){
			string command = "rm -rf " + this->_outputDirTemp + "/solid/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/temp/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/count_synchro/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/merge_synchro/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/stats/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/job_count/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/job_merge/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/kmercount_per_partition/";
			system(command.c_str());
			command = "rm -rf " + this->_outputDirTemp + "/input/";
			system(command.c_str());


			command = "rm " + this->_outputDirTemp + "/config.h5";
			system(command.c_str());
			command = "rm " + this->_outputDirTemp + "/datasetIds";
			system(command.c_str());
			//cout << command << endl;
			//System::file().rmdir(this->_outputDirTemp);

			//System::file().mkdir(this->_outputDirTemp + "/solid/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/temp/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/log/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/count_synchro/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/merge_synchro/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/stats/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/job_count/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/job_merge/", -1);
			//System::file().mkdir(this->_outputDirTemp + "/kmercount_per_partition/", -1);
		}
	}

	void parseArgs() {

		SimkaAlgorithm<span>::parseArgs();

		if(this->_options->get(STR_SIMKA_JOB_COUNT_FILENAME) || this->_options->get(STR_SIMKA_JOB_MERGE_FILENAME) || this->_options->get(STR_SIMKA_JOB_COUNT_COMMAND) || this->_options->get(STR_SIMKA_JOB_MERGE_COMMAND)){
			_isClusterMode = true;
			_jobCountFilename = this->_options->getStr(STR_SIMKA_JOB_COUNT_FILENAME);
			_jobMergeFilename = this->_options->getStr(STR_SIMKA_JOB_MERGE_FILENAME);
			_jobCountCommand = this->_options->getStr(STR_SIMKA_JOB_COUNT_COMMAND);
			_jobMergeCommand = this->_options->getStr(STR_SIMKA_JOB_MERGE_COMMAND);

			if(! this->_options->get(STR_SIMKA_NB_JOB_COUNT) || this->_options->get(STR_SIMKA_NB_JOB_MERGE)){
				//cout << endl;
				cout << "Cluster mode enable. Be sure to set correctly the following arguments if you have any job submission constraints:" << endl;
				cout << "\t" << STR_SIMKA_NB_JOB_COUNT << " : the maximum number of simultaneous couting" << endl; //job (each job will use up to " << STR_NB_CORES << " cores and " << STR_MAX_MEMORY << " MB memory)" << endl;
				cout << "\t" << STR_SIMKA_NB_JOB_MERGE << " : the maximum number of simultaneous merging job" << endl; // (each job will use up to 1 core and " << STR_MAX_MEMORY << " MB memory)" << endl;
				//cout << endl;
			}




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


		}
		else{
			_isClusterMode = false;
		}

        if (this->_options->getInt(STR_MAX_MEMORY) < 2000)
        {
            std::cout << "WARNING: running Simka with low memory is risky. Simka may hang because of that. Consider running with -max-memory X where X > 2000" << std::endl;
        }



		if (this->_options->getInt(STR_MAX_MEMORY) < 500)
        {
            std::cout << "Please run Simka with higher memory usage than 500 MB" << std::endl;
            exit(1);
        }


		//_isClusterMode = true;

		//if(this->_options->get(STR_SIMKA_CLUSTER_MODE)){

			//cout << "cluster mode activated" << endl;
			//cout << "\t-max-memory = memory per job" << endl;
			//cout << "\t-nb-cores = cores per job" << endl;
			//cout << endl;


			//if(_isClusterMode)
			//	_maxJobMerge = max((int)_maxJobMerge, (int)30);
			//else{
			//	_maxJobMerge = max((int)_maxJobMerge, (int)30);
			//}



		//}


		/*
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
		*/
	}


	void setup(){
		SimkaAlgorithm<span>::setup();

		createDirs();
		layoutInputFilename();

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
		//System::file().mkdir(this->_outputDirTemp + "/solid/merged/", -1);
		System::file().mkdir(this->_outputDirTemp + "/temp/", -1);
		System::file().mkdir(this->_outputDirTemp + "/log/", -1);
		System::file().mkdir(this->_outputDirTemp + "/count_synchro/", -1);
		System::file().mkdir(this->_outputDirTemp + "/merge_synchro/", -1);
		System::file().mkdir(this->_outputDirTemp + "/stats/", -1);
		System::file().mkdir(this->_outputDirTemp + "/job_count/", -1);
		System::file().mkdir(this->_outputDirTemp + "/job_merge/", -1);
		System::file().mkdir(this->_outputDirTemp + "/kmercount_per_partition/", -1);

	}

	void createConfig(){

		size_t maxCores = this->_nbCores;
		size_t maxMemory = this->_maxMemory;
		size_t minMemoryPerJobMB = 500;


		if(this->_options->get(STR_SIMKA_NB_JOB_COUNT)){
			_maxJobCount = this->_options->getInt(STR_SIMKA_NB_JOB_COUNT);
			//maxCores = _maxJobCount; //TO REMOVE WHEN BUG IN DISPATCHER IS RESOLVED
		}
		else{
			size_t maxjob_byCore = min(maxCores/2, this->_nbBanks);
			//size_t maxjob_byCore = min(maxCores, this->_nbBanks);  //TO REMOVE WHEN BUG IN DISPATCHER IS RESOLVED

			maxjob_byCore = max(maxjob_byCore, (size_t)1);

			size_t maxjob_byMemory = maxMemory/minMemoryPerJobMB;
			maxjob_byMemory = max(maxjob_byMemory, (size_t) 1);

			size_t maxJobs = min(maxjob_byCore, maxjob_byMemory);
			_maxJobCount = maxJobs;
			//maxCores = _maxJobCount; //TO REMOVE WHEN BUG IN DISPATCHER IS RESOLVED

		}
		//_maxJobCount = 1;

		if(this->_options->get(STR_SIMKA_NB_JOB_MERGE)){
			_maxJobMerge = this->_options->getInt(STR_SIMKA_NB_JOB_MERGE);
		}
		else{
			_maxJobMerge = maxCores;
			/*
			if(this->_computeComplexDistances && this->_computeSimpleDistances){
				_maxJobMerge = max((size_t)maxCores/4, (size_t)1);
			}
			else if(this->_computeSimpleDistances){
				_maxJobMerge = max((size_t)maxCores/2, (size_t)1);
			}
			else if(this->_computeComplexDistances){
				_maxJobMerge = max((size_t)maxCores/3, (size_t)1);
			}
			else{
				_maxJobMerge = maxCores;
			}*/
		}

		_maxJobCount = min(_maxJobCount, maxCores);
		_maxJobMerge = min(_maxJobMerge, maxCores);

		_coresPerJob = maxCores / _maxJobCount;
		_coresPerJob = max((size_t)1, _coresPerJob);

		_memoryPerJob = maxMemory / _maxJobCount;
		_memoryPerJob = max(_memoryPerJob, (size_t)minMemoryPerJobMB);

		_coresPerMergeJob = maxCores / _maxJobMerge;
		_coresPerMergeJob = max((size_t)1, _coresPerMergeJob);

		cout << endl;
		cout << "Maximum ressources used by Simka: " << endl;
		cout << "\t - " << _maxJobCount << " simultaneous processes for counting the kmers (per job: " << _coresPerJob << " cores, " << _memoryPerJob << " MB memory)" << endl;
		cout << "\t - " << _maxJobMerge << " simultaneous processes for merging the kmer counts (per job: " << _coresPerMergeJob << " cores, memory undefined)" << endl;
		cout << endl;



		//_coresPerJob = this->_nbCores;
		//_memoryPerJob = this->_maxMemory;















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

		this->_options->setInt(STR_MAX_MEMORY, _memoryPerJob - _memoryPerJob/3);
		this->_options->setInt(STR_NB_CORES, _coresPerJob);

	    Storage* storage = 0;
        storage = StorageFactory(STORAGE_HDF5).create (filename, true, false);
        LOCAL (storage);





        size_t chosenBankId = 0;
    	SimkaSequenceFilter dummyFilter(0, 0);
    	//vector<SimkaBankFiltered<SimkaSequenceFilter>*> banksToDelete;

    	string inputDir = this->_outputDirTemp + "/input/";
    	u_int64_t maxPart = 0;
    	for (size_t i=0; i<this->_nbBanks; i++){

    		IBank* bank = Bank::open(inputDir + this->_bankNames[i]);
    		LOCAL(bank);

            // Check that data files are not empty
            int64_t nbItems = bank->estimateNbItems();
            if (nbItems==0) {
                cerr << "ERROR: Dataset is empty: " << this->_bankNames[i] << endl;
                exit(1);
            }

    		//size_t nbBank_ = bank->getCompositionNb();
    		SimkaBankTemp* simkaBank = new SimkaBankTemp(bank, this->_maxNbReads*this->_nbBankPerDataset[i]);
    		//banksToDelete.push_back(simkaBank);
    		ConfigurationAlgorithm<span> testConfig(simkaBank, this->_options);
    		testConfig.execute();

    		size_t part = testConfig.getConfiguration()._nb_partitions;
    		if(part > maxPart){
    			maxPart = part;
    			chosenBankId = i;
    		}

    		//delete simkaBank;
/*
    		u_int64_t nbReads = bank->estimateNbItems();
    		nbReads /= _nbBankPerDataset[i];
    		totalReads += nbReads;
    		if(nbReads < minReads){
    			minReads = nbReads;
    			//_smallerBankId = _bankNames[i];
    		}
    		if(nbReads > maxReads){
    			maxReads = nbReads;
    			_largerBankId = _bankNames[i];
    		}*/

    	}


    	//maxPart += 2;
    	//for(size_t i=0; i<banksToDelete.size(); i++)
    	//	delete banksToDelete[i];



		this->_options->setInt(STR_MAX_MEMORY, _memoryPerJob);

    	IBank* inputbank = Bank::open(this->_banksInputFilename);
    	LOCAL(inputbank);

		IBank* bank = Bank::open(this->_outputDirTemp + "/input/" + this->_bankNames[chosenBankId]);
		LOCAL(bank);

		//IBank* bank = Bank::open(_outputDirTemp + "/input/" + _bankNames[0]);
		//LOCAL(inputbank);

        //bank->finalize();

		//u_int64_t nbSeqs = 1;
        //IBank* sampleBank = new SimkaBankSample(bank, nbSeqs);
		//SortingCountAlgorithm<span> sortingCount (sampleBank, this->_options);

		//SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();

		//sortingCount.addProcessor (proc);

		// We launch the algorithm
		//sortingCount.execute();

		//Configuration config = sortingCount.getConfig();

		ConfigurationAlgorithm<span> testConfig1(inputbank, this->_options);
		testConfig1.execute();
		Configuration config1 = testConfig1.getConfiguration();


		ConfigurationAlgorithm<span> testConfig2(bank, this->_options);
		testConfig2.execute();
		Configuration config2 = testConfig2.getConfiguration();


        //IBank* inputbank = Bank::open(this->_banksInputFilename);
		//LOCAL(inputbank);
		//IBank* sampleBank2 = new SimkaBankSample(inputbank, nbSeqs);
		//SortingCountAlgorithm<span> sortingCount2 (sampleBank2, this->_options);
		//SimkaNullProcessor<span>* proc2 = new SimkaNullProcessor<span>();
		//sortingCount2.addProcessor (proc2);
		//sortingCount2.execute();
		//Configuration config2 = sortingCount2.getConfig();
		//cout << config2._nb_partitions << endl;


		_nbPartitions = max((size_t)maxPart, (size_t)_maxJobMerge);
		//_nbPartitions = max(_nbPartitions, (size_t)32);

		cout << "Nb partitions: " << _nbPartitions << " partitions" << endl << endl << endl;
		//_nbPartitions = max((int)_nbPartitions, (int)30);

		config1._nb_partitions = _nbPartitions;
		config2._nb_partitions = _nbPartitions;

        RepartitorAlgorithm<span> repart (inputbank, storage->getGroup(""), config1);
        repart.execute ();


		uint64_t memoryUsageCachedItems;
		config2._nb_cached_items_per_core_per_part = 1 << 8; // cache at least 256 items (128 here, then * 2 in the next while loop)
		do
		{
			config2._nb_cached_items_per_core_per_part *= 2;
			memoryUsageCachedItems = 1LL * config2._nb_cached_items_per_core_per_part *config2._nb_partitions * config2._nbCores * sizeof(Type);
		}
		while (memoryUsageCachedItems < config2._max_memory * MBYTE / 10);
		/*
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
		}*/

		//IBank* inputbank = Bank::open(this->_banksInputFilename);
		//LOCAL(inputbank);
		//ConfigurationAlgorithm<span> inputconfig(inputbank, this->_options);
		//inputconfig.execute();
        //RepartitorAlgorithm<span> repart (inputbank, storage->getGroup(""), config);
        //repart.execute ();
        //setRepartitor (new Repartitor(storage->getGroup("minimizers")));
		//SortingCountAlgorithm<span> sortingCount (sampleBank, _options);




		config2.save(storage->getGroup(""));
		//sortingCount.getRepartitor()->save(storage->getGroup(""));
		//delete sampleBank;

	    //setStorage (storage);

		//delete storage;
		//sampleBank->forget();
	}

	void removeMergeSynchro(){

	    for (size_t i=0; i<this->_bankNames.size(); i++){
			string finishFilename = this->_outputDirTemp + "/merge_synchro/" +  this->_bankNames[i] + ".ok";
			if(System::file().doesExist(finishFilename)) System::file().remove(finishFilename);
	    }
	}

	void printCountInfo(){

		char * pEnd;
		vector<u_int64_t> kmerPerParts(_nbPartitions, 0);


		for(size_t i=0; i<this->_bankNames.size(); i++){
			//cout << filename << endl;

			string line;
			size_t currentPart = 0;
			ifstream file((this->_outputDirTemp + "/kmercount_per_partition/" +  this->_bankNames[i] + ".txt").c_str());
			size_t j = 0;
			while(getline(file, line)){
				if(line == "") continue;
				kmerPerParts[j] += strtoull(line.c_str(), NULL, 10);
				j += 1;
			}
			file.close();
    	}

		cout << endl << endl << "Kmer repartition" << endl;
		for(size_t i=0; i<kmerPerParts.size(); i++){
			cout <<  "\t" << i << ":\t" << kmerPerParts[i] << endl;
		}
		cout << endl << endl;
	}

	void count(){

		cout << endl << "Counting k-mers... (log files are " + this->_outputDirTemp + "/log/count_*)" << endl;

	    for (size_t i=0; i<_nbPartitions; i++){
	    	//System::file().mkdir(this->_outputDirTemp + "/solid/merged/part_" + Stringify::format("%i", i), -1);
	    	System::file().mkdir(this->_outputDirTemp + "/solid/part_" + Stringify::format("%i", i), -1);
	    }

		vector<string> commands;

		_progress = new ProgressSynchro (
			this->createIteratorListener (this->_bankNames.size(), "Counting datasets"),
			System::thread().newSynchronizer());
		_progress->init ();

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<this->_bankNames.size(); i++){

			string logFilename = this->_outputDirTemp + "/log/count_" + this->_bankNames[i] + ".txt";

			string finishFilename = this->_outputDirTemp + "/count_synchro/" +  this->_bankNames[i] + ".ok";
			if(System::file().doesExist(finishFilename)){
				_progress->inc(1);
				cout << "\t" << this->_bankNames[i] << " already counted (remove file " << finishFilename << " to count again)" << endl;
				continue;
			}
			//else{

			string tempDir = this->_outputDirTemp + "/temp/" + this->_bankNames[i];

			string command = "nohup " + _execDir + "simkaCountProcess " + _execDir + "simkaCount ";
			command += " " + string(STR_KMER_SIZE) + " " + SimkaAlgorithm<>::toString(this->_kmerSize);
			command += " " + string("-out-tmp-simka") + " " + this->_outputDirTemp;
			command += " " + string("-out-tmp") + " " + tempDir;
			command += " -bank-name " + this->_bankNames[i];
			command += " -bank-index " + SimkaAlgorithm<>::toString(i);
			command += " -nb-datasets " + SimkaAlgorithm<>::toString(this->_nbBankPerDataset[i]);
			command += " " + string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(_memoryPerJob);
			command += " " + string(STR_NB_CORES) + " " + SimkaAlgorithm<>::toString(_coresPerJob);
			command += " " + string(STR_URI_INPUT) + " dummy ";
			command += " " + string(STR_KMER_ABUNDANCE_MIN) + " " + SimkaAlgorithm<>::toString(this->_abundanceThreshold.first);
			command += " " + string(STR_KMER_ABUNDANCE_MAX) + " " + SimkaAlgorithm<>::toString(this->_abundanceThreshold.second);
			command += " " + string(STR_SIMKA_MIN_READ_SIZE) + " " + SimkaAlgorithm<>::toString(this->_minReadSize);
			command += " " + string(STR_SIMKA_MIN_READ_SHANNON_INDEX) + " " + Stringify::format("%f", this->_minReadShannonIndex);
			command += " " + string(STR_SIMKA_MAX_READS) + " " + SimkaAlgorithm<>::toString(this->_maxNbReads);
			command += " -nb-partitions " + SimkaAlgorithm<>::toString(_nbPartitions);
			//command += " -verbose " + Stringify::format("%d", this->_options->getInt(STR_VERBOSE));
			command += " >> " + logFilename + " 2>&1";

			filenameQueue.push_back(this->_bankNames[i]);
			System::file().mkdir(tempDir, -1);

			string str = "Counting dataset " + SimkaAlgorithm<>::toString(i) + "\n";
			str += "\t" + command + "\n\n\n";
			system(("echo \"" + str + "\" > " + logFilename).c_str());

			//cout << "Counting dataset " << i << endl;
			//cout << "\t" << command << endl;

			removeMergeSynchro();



			//_progress->inc(1);
			//nanosleep((const struct timespec[]){{0, 10000000L}}, NULL);


			if(_isClusterMode){
				string jobFilename = this->_outputDirTemp + "/job_count/job_count_" + SimkaAlgorithm<>::toString(i) + ".bash";
				IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
				system(("chmod 755 " + jobFilename).c_str());
				string jobCommand = _jobCountContents + '\n' + '\n';
				jobCommand += command;


				jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
				jobFile->flush();
				string submitCommand = _jobCountCommand + " " + jobFile->getPath();
				delete jobFile;
				system(submitCommand.c_str());
			}
			else{
				command += " &";
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
						//sleep(1);
						const struct timespec t[] {{ 0, 1000000000L}};
                        nanosleep(t, NULL);
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





		cout << endl << "Merging k-mer counts and computing distances... (log files are " + this->_outputDirTemp + "/log/merge_*)" << endl;

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

			string logFilename = this->_outputDirTemp + "/log/merge_" + datasetId + ".txt";

			if(System::file().doesExist(finishFilename)){
				_progress->inc(1);
				cout << "\t" << datasetId << " already merged (remove file " << finishFilename << " to merge again)" << endl;
			}
			else{
				//if(System::file().doesExist(finishFilename)){
				//	System::file().remove(finishFilename);
				//	cout << "\t" << _bankNames[i] << " already  (remove file " << finishFilename << " to count again)" << endl;
				//}
				//else{

				filenameQueue.push_back(datasetId);

				string command = "nohup " + _execDir + "simkaMerge ";
				command += " " + string(STR_KMER_SIZE) + " " + SimkaAlgorithm<>::toString(this->_kmerSize);
				command += " " + string(STR_URI_INPUT) + " " + this->_inputFilename;
				command += " " + string("-out-tmp-simka") + " " + this->_outputDirTemp;
				command += " -partition-id " + SimkaAlgorithm<>::toString(i);
				command += " " + string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(this->_maxMemory / this->_nbCores);
				command += " " + string(STR_NB_CORES) + " " + SimkaAlgorithm<>::toString(_coresPerMergeJob);
				command += " " + string(STR_SIMKA_MIN_KMER_SHANNON_INDEX) + " " + Stringify::format("%f", this->_minKmerShannonIndex);
				command += " -verbose " + Stringify::format("%d", this->_options->getInt(STR_VERBOSE));
				if(this->_computeSimpleDistances) command += " " + string(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
				if(this->_computeComplexDistances) command += " " + string(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);
				command += " >> " + logFilename + " 2>&1";
				//SimkaDistanceParam distanceParams(this->_options);
				//if(distanceParams._computeBrayCurtis) command += " " + STR_SIMKA_DISTANCE_BRAYCURTIS + " ";
				//if(distanceParams._computeCanberra) command += " " + STR_SIMKA_DISTANCE_CANBERRA + " ";
				//if(distanceParams._computeChord) command += " " + STR_SIMKA_DISTANCE_CHORD + " ";
				//if(distanceParams._computeHellinger) command += " " + STR_SIMKA_DISTANCE_HELLINGER + " ";
				//if(distanceParams._computeKulczynski) command += " " + STR_SIMKA_DISTANCE_KULCZYNSKI + " ";


				string str = "Merging partition " + SimkaAlgorithm<>::toString(i) + "\n";
				str += "\t" + command + "\n\n\n";
				system(("echo \"" + str + "\" > " + logFilename).c_str());


				if(_isClusterMode){
					string jobFilename = this->_outputDirTemp + "/job_merge/job_merge_" + SimkaAlgorithm<>::toString(i) + ".bash";
					IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
					system(("chmod 755 " + jobFilename).c_str());
					string jobCommand = _jobMergeContents + '\n' + '\n';
					jobCommand += command;


					jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
					jobFile->flush();
					string submitCommand = _jobMergeCommand + " " + jobFile->getPath();
					delete jobFile;
					system(submitCommand.c_str());
				}
				else{
					command += " &";
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

	    //cout << nbJobs << endl;

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


	    //cout << nbJobs << endl;

	    _progress->finish();
	    delete _progress;
	}

	/*
	void getCountInfo(SimkaStatistics& mainStats){

    	for(size_t i=0; i<this->_nbBanks; i++){
    		string name = this->_bankNames[i];
    		string countFilename = this->_outputDirTemp + "/count_synchro/" +  name + ".ok";

    		string line;
	    	ifstream file(countFilename.c_str());
	    	vector<string> lines;
			while(getline(file, line)){
				if(line == "") continue;
				lines.push_back(line);
			}
			file.close();

			u_int64_t nbReads = strtoull(lines[0].c_str(), NULL, 10) ;
			mainStats._datasetNbReads[i] = nbReads;
			mainStats._nbSolidDistinctKmersPerBank[i] = strtoull(lines[1].c_str(), NULL, 10);
			//cout << mainStats._nbSolidDistinctKmersPerBank[i] << endl;
			mainStats._nbSolidKmersPerBank[i] = strtoull(lines[2].c_str(), NULL, 10);
			mainStats._chord_sqrt_N2[i] = sqrt(strtoull(lines[3].c_str(), NULL, 10));
    	}

	}*/

	void stats(){
		cout << endl << "Computing stats..." << endl;
		//cout << this->_nbBanks << endl;

		//u_int64_t nbKmers = 0;

		//SimkaDistanceParam distanceParams(this->_options);
		SimkaStatistics mainStats(this->_nbBanks, this->_computeSimpleDistances, this->_computeComplexDistances, this->_outputDirTemp, this->_bankNames);

		for(size_t i=0; i<_nbPartitions; i++){

			string filename = this->_outputDirTemp + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".gz";
			//Storage* storage = StorageFactory(STORAGE_HDF5).load (this->_outputDirTemp + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".stats");
			//LOCAL (storage);

			SimkaStatistics stats(this->_nbBanks, this->_computeSimpleDistances, this->_computeComplexDistances, this->_outputDirTemp, this->_bankNames);
			stats.load(filename);

			//cout << stats._nbDistinctKmers << "   " << stats._nbKmers << endl;
			mainStats += stats;

			//nbKmers += stats._nbKmers;
		}

		//cout << "Nb kmers: " << nbKmers << endl;

		//getCountInfo(mainStats);
		//for(size_t i=0; i<this->_nbBanks; i++){
		//	cout << mainStats._nbSolidDistinctKmersPerBank[i] << endl;
		//}
		mainStats.outputMatrix(this->_outputDir, this->_bankNames);

#//ifdef PRINT_STATS
		if(this->_options->getInt(STR_VERBOSE) != 0) mainStats.print();
#//endif
	}


	//u_int64_t _maxMemory;
	//size_t _nbCores;
	size_t _memoryPerJob;
	size_t _coresPerJob;
	size_t _coresPerMergeJob;

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

    string _execDir;
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

	SimkaPotara(const string& execFilename);
	void execute();

	string _execFilename;
};

#endif
