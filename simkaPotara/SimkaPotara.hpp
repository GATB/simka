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

//#define CLUSTER
//#define SERIAL
#define SLEEP_TIME_SEC 1


const string STR_SIMKA_JOB_COUNT_COMMAND = "-job-count-cmd";
const string STR_SIMKA_JOB_MERGE_COMMAND = "-job-merge-cmd";
const string STR_SIMKA_JOB_COUNT_FILENAME = "-job-count-file";
const string STR_SIMKA_JOB_MERGE_FILENAME = "-job-merge-file";
const string STR_SIMKA_NB_JOB_COUNT = "-max-job-count";
const string STR_SIMKA_NB_JOB_MERGE = "-max-job-merge";
const string STR_SIMKA_NB_PARTITIONS = "-nb-partitions";

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
class SimkaPotaraAlgorithm : public Algorithm{
public:


	SimkaPotaraAlgorithm(IProperties* options):
	Algorithm("SimkaPotara", 1, options)
	{

		_options = options;

		_inputFilename = _options->getStr(STR_URI_INPUT);
		_outputDir = _options->getStr(STR_URI_OUTPUT);
		_outputDirTemp = _options->getStr(STR_URI_OUTPUT_TMP);

		_maxNbReads = _options->getInt(STR_SIMKA_MAX_READS);
		//_maxJobCount = _options->getInt(STR_SIMKA_NB_JOB_COUNT);
		//_maxJobMerge = _options->getInt(STR_SIMKA_NB_JOB_MERGE);
		//_jobCountFilename = _options->getStr(STR_SIMKA_JOB_COUNT_FILENAME);
		//_jobMergeFilename = _options->getStr(STR_SIMKA_JOB_MERGE_FILENAME);
		//_jobCountCommand = _options->getStr(STR_SIMKA_JOB_COUNT_COMMAND);
		//_jobMergeCommand = _options->getStr(STR_SIMKA_JOB_MERGE_COMMAND);

		_nbAskedPartitions = _options->getInt(STR_SIMKA_NB_PARTITIONS);

		_maxMemory = _options->getInt(STR_MAX_MEMORY);
		_nbCores = _options->getInt(STR_NB_CORES);

		_maxJobMerge = _nbCores;

		size_t maxCoreCount = (2*System::info().getNbCores()) / 3;
		size_t nbCoresCount = min(maxCoreCount, _nbCores);

		u_int64_t minMemory = 2000;
		size_t maxJobCountTemp = _maxMemory/minMemory;
		_maxJobCount = min(nbCoresCount, maxJobCountTemp);
		_memoryPerJob = _maxMemory / _maxJobCount;
		_coresPerJob = ceil(nbCoresCount / (float)_maxJobCount);


		cout << "Nb jobs in parallel: " << _maxJobCount << endl;
		cout << "Cores per jobs: " << _coresPerJob << endl;
		cout << "Memory per jobs: " << _memoryPerJob << endl;
		//string solidFilename = _outputDir + "/solid/" +  p.bankName + suffix + ".h5";

		//cout << "SimkaFusion constructor       " << _outputDirTempFilter << endl;




	}

	~SimkaPotaraAlgorithm(){

	}


	void execute(){

		if(!System::file().doesExist(_outputDir)){
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
		        std::cout << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
		        return;
			}
		}

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
		ifstream inputFile(_inputFilename.c_str());
		IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		string line;
		string linePart;
		vector<string> linePartList;

		string bankFileContents = "";

		u_int64_t lineIndex = 0;

		while(getline(inputFile, line)){

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

				bankFileContents += subBankFilename + "\n";
				_nbBankPerDataset.push_back(linePartList.size() - 1); //linePartList.size() - 1 = nb sub banks
				//_nbReadsPerDataset.push_back(ceil(_maxNbReads / (float)()));
			}

			lineIndex += 1;
		}

		bankFileContents.erase(bankFileContents.size()-1);
		//bankFileContents.pop_back(); // "remove last /n

		bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);

		inputFile.close();
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

		ifstream inputFile(_inputFilename.c_str());


		//Hold dataset ids
		string datasetIdFilename = _outputDirTempFilter + "/" + "datasetIds";
		IFile* datasetIdFile = System::file().newFile(datasetIdFilename, "wb");


		string line;
		string linePart;
		vector<string> linePartList;

		//string bankFileContents = "";

		u_int64_t lineIndex = 0;

		while(getline(inputFile, line)){

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
		inputFile.close();
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
			if(_options->getInt(STR_VERBOSE) != 0)
				cout << "-maxNbReads is not defined. Simka will estimating it..." << endl;
			//_maxNbReads = bank->estimateNbItems() / _nbBanks;
			//_maxNbReads -= (_maxNbReads/10);
			u_int64_t minReads = -1;
			for (size_t i=0; i<_nbBanks; i++){
				u_int64_t nbReads = bank->estimateNbItemsBanki(i);
				if(nbReads < minReads) minReads = nbReads;
			}
			_maxNbReads = minReads;
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

		string suffix = "";
		suffix += "m" + _options->getStr(STR_SIMKA_MIN_READ_SIZE);
		suffix += "_s" + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX);
		suffix += "_n" + SimkaAlgorithm<>::toString(_maxNbReads);
		suffix += "_p" + SimkaAlgorithm<>::toString(_nbAskedPartitions);
		_outputDirTempFilter = _outputDirTemp + "/" + suffix + "/";

		System::file().mkdir(_outputDirTemp, -1);
		System::file().mkdir(_outputDirTempFilter, -1);
		System::file().mkdir(_outputDirTempFilter + "/input/", -1);
		System::file().mkdir(_outputDirTempFilter + "/solid/", -1);
		System::file().mkdir(_outputDirTempFilter + "/temp/", -1);
		System::file().mkdir(_outputDirTempFilter + "/count_synchro/", -1);
		System::file().mkdir(_outputDirTempFilter + "/merge_synchro/", -1);
		System::file().mkdir(_outputDirTempFilter + "/stats/", -1);
		System::file().mkdir(_outputDirTempFilter + "/job_count/", -1);
		System::file().mkdir(_outputDirTempFilter + "/job_merge/", -1);

#ifdef CLUSTER
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
#endif
	}

	void createConfig(){

		string filename = _outputDirTempFilter + "/" + "config.h5";
		if(System::file().doesExist(filename)){
			cout << "\t" << " config already exists (remove file " << filename << " to config again)" << endl;


			Storage* storage = StorageFactory(STORAGE_HDF5).load (filename);
			LOCAL (storage);
			Configuration* config = new Configuration();
			config->load(storage->getGroup(""));
			_nbPartitions = config->_nb_partitions;

			return;
		}

		_options->setInt(STR_MAX_MEMORY, _memoryPerJob);
		_options->setInt(STR_NB_CORES, _coresPerJob);

	    Storage* storage = 0;
        storage = StorageFactory(STORAGE_HDF5).create (filename, true, false);
        LOCAL (storage);

        IBank* bank = Bank::open(_outputDirTempFilter + "/input/" + _bankNames[0]);
        bank->finalize();
        IBank* sampleBank = new SimkaBankSample(bank, _maxNbReads/3);
		SortingCountAlgorithm<span> sortingCount (sampleBank, _options);

		SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();

		sortingCount.addProcessor (proc);

		// We launch the algorithm
		sortingCount.execute();

		Configuration config = sortingCount.getConfig();
		if(_nbAskedPartitions == 0){
			_nbPartitions = config._nb_partitions;
		}
		else{
			_nbPartitions = _nbAskedPartitions;
			config._nb_partitions = _nbPartitions;
		}


        RepartitorAlgorithm<span> repart (bank, storage->getGroup(""), config);
        repart.execute ();
        //setRepartitor (new Repartitor(storage->getGroup("minimizers")));
		//SortingCountAlgorithm<span> sortingCount (sampleBank, _options);




		config.save(storage->getGroup(""));
		//sortingCount.getRepartitor()->save(storage->getGroup(""));


	    //setStorage (storage);

		//delete storage;
		//sampleBank->forget();
	}


	void count(){

		vector<vector<string> > commands;

		_progress = new ProgressSynchro (
			createIteratorListener (_bankNames.size(), "Counting datasets"),
			System::thread().newSynchronizer());
		_progress->init ();

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;


			/*
		id_t pids[_maxJobCount];
		//int pipe_fd[_maxJobCount];
		for(int child=0;child<_maxJobCount;child++) {
			cout << child << endl;
			 //int pipe[2];
			 int ret;
			 ret = fork();
			 if(ret) {
			  pids[child] = ret;
			  cout << "ha:   " << ret << endl;
			 } else {
				 cout << "child" << endl;
			  //exec(...);
			 }
		}
*/

	    for (size_t i=0; i<_bankNames.size(); i++){

	    	cout << endl << endl << "\tDataset " << i << endl;

			string finishFilename = _outputDirTempFilter + "/count_synchro/" +  _bankNames[i] + ".ok";
			if(System::file().doesExist(finishFilename)){
				cout << "\t" << _bankNames[i] << " already counted (remove file " << finishFilename << " to count again)" << endl;
			}
			else{

				string tempDir = _outputDirTempFilter + "/temp/" + _bankNames[i];

				vector<string> command;
				command.push_back(_bankNames[i]);
				command.push_back(tempDir);
				/*
				command.push_back();
				command.push_back(+ " " + );

				command.push_back(+ " " + );
				command.push_back( + );
				command.push_back(string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(_memoryPerJob));
				command.push_back(string(STR_NB_CORES) + " " + SimkaAlgorithm<>::toString(_coresPerJob));
				command.push_back(string() + " dummy ");
				command.push_back(string(STR_KMER_ABUNDANCE_MIN) + " " + _options->getStr(STR_KMER_ABUNDANCE_MIN));
				command.push_back(string(STR_SIMKA_MIN_READ_SIZE) + " " + _options->getStr(STR_SIMKA_MIN_READ_SIZE));
				command.push_back(string(STR_SIMKA_MIN_READ_SHANNON_INDEX) + " " + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX));
				command.push_back(string(STR_SIMKA_MAX_READS) + " " + SimkaAlgorithm<>::toString(_nbReadsPerDataset[i]));
				*/
				commands.push_back(command);
			}
	    }

	    vector<pid_t> pids(commands.size());
	    size_t jobDone = 0;
	    vector<bool> isJobFinished(commands.size());
	    //size_t nbJobs = 0;

	    for(size_t i=0; i<commands.size(); i++){
			 int ret;
			ret = fork();
			if(ret) {
				pids[i] = ret;
				nbJobs += 1;
				cout << "-------------:   " << ret << "   " << i << "   " << nbJobs << endl;

				if(nbJobs >= _maxJobCount){

					//while(true){
					//bool ok = false;
					for(size_t j=0; j<=i; j++){

						if(isJobFinished[j]) continue;

						int status;
						waitpid(pids[j], &status, 0);
						cout << "job finiashed " << j << endl;

						isJobFinished[j] = true;
						nbJobs -= 1;
						//jobDone += 1;
						//ok = true;
					}
						//if(ok) break;
					//}
				}

				cout << nbJobs << endl;
				//if(i == commands.size()-1)
				//    waitpid(pids[i], &status, 0);
				//nbJobs -= 1;
				//}

			}
			else {
				vector<string> command = commands[i];
				cout << "child" << endl;
				execl(
						"./simkaCount",
						(string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE)).c_str(),
						string("-out-tmp-simka") .c_str(),
						_outputDirTempFilter.c_str(),
						string("-out-tmp") .c_str(),
						command[1].c_str(),
						"-bank-name" ,
						command[0].c_str(),
						STR_URI_INPUT,
						"dummy",
						STR_SIMKA_MAX_READS.c_str(),
						SimkaAlgorithm<>::toString(_maxNbReads).c_str(),
						STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(),
						"0",
						STR_SIMKA_MIN_READ_SIZE.c_str(),
						"0",
						STR_KMER_ABUNDANCE_MIN,
						_options->getStr(STR_KMER_ABUNDANCE_MIN).c_str(),
						"-verbose",
						"0",
						NULL
				);
			}
	    }


		for(size_t j=0; j<commands.size(); j++){
			int status;
			waitpid(pids[j], &status, 0);
			//ok = true;
		}

	    /*
				filenameQueue.push_back(_bankNames[i]);


				string tempDir = _outputDirTempFilter + "/temp/" + _bankNames[i];

				string command = "";
#ifndef CLUSTER
				//command += "nohup";
#endif
				command += "./simkaCount ";
				command += " " + string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE);
				command += " " + string("-out-tmp-simka") + " " + _outputDirTempFilter;
				command += " " + string("-out-tmp") + " " + tempDir;
				command += " -bank-name " + _bankNames[i];
				command += " " + string(STR_MAX_MEMORY) + " " + SimkaAlgorithm<>::toString(_memoryPerJob);
				command += " " + string(STR_NB_CORES) + " " + SimkaAlgorithm<>::toString(_coresPerJob);
				command += " " + string(STR_URI_INPUT) + " dummy ";
				command += " " + string(STR_KMER_ABUNDANCE_MIN) + " " + _options->getStr(STR_KMER_ABUNDANCE_MIN);
				command += " " + string(STR_SIMKA_MIN_READ_SIZE) + " " + _options->getStr(STR_SIMKA_MIN_READ_SIZE);
				command += " " + string(STR_SIMKA_MIN_READ_SHANNON_INDEX) + " " + _options->getStr(STR_SIMKA_MIN_READ_SHANNON_INDEX);
				command += " " + string(STR_SIMKA_MAX_READS) + " " + SimkaAlgorithm<>::toString(_nbReadsPerDataset[i]);
#ifndef CLUSTER
				command += " &";
#endif

				//command = _cmdJobCount + " " + command;

#ifdef CLUSTER
				string jobFilename = _outputDirTempFilter + "/job_count/job_count_" + _bankNames[i] + ".bash";
				IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
				string jobCommand = _jobCountContents + '\n' + '\n';
				jobCommand += command;

				cout << "\t" << jobCommand << endl;

				jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
				jobFile->flush();
				string submitCommand = _jobCountCommand + " " + jobFile->getPath();
				delete jobFile;
				system(submitCommand.c_str());
#else
				cout << "\t" << command << endl;
				//execl("./simkaCount",
				//		(string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE)).c_str(),
				//	NULL);
				system(command.c_str());
#endif
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


	    cout << nbJobs << endl;*/

	    _progress->finish();
	    delete _progress;
	}

	void waitForJob(){

	}

	void merge(){

		//cout << "lala" << endl;
		_progress = new ProgressSynchro (
			createIteratorListener (_nbPartitions, "Merging datasets"),
			System::thread().newSynchronizer());
		_progress->init ();
		//cout << "lala" << endl;

		vector<string> filenameQueue;
		vector<string> filenameQueueToRemove;
		size_t nbJobs = 0;

	    for (size_t i=0; i<_nbPartitions; i++){

	    	cout << endl << endl <<  "\tDataset " << i << endl;

	    	string datasetId = SimkaAlgorithm<>::toString(i);
			string finishFilename = _outputDirTempFilter + "/merge_synchro/" +  datasetId + ".ok";

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

				string command = "";

	#ifndef CLUSTER
				//command += "nohup";
	#endif
				command += "./simkaMerge ";
				command += " " + string(STR_KMER_SIZE) + " " + _options->getStr(STR_KMER_SIZE);
				command += " " + string(STR_URI_INPUT) + " " + _inputFilename;
				command += " " + string("-out-tmp-simka") + " " + _outputDirTempFilter;
				command += " -partition-id " + SimkaAlgorithm<>::toString(i);
				command += " " + string(STR_MAX_MEMORY) + " " + _options->getStr(STR_MAX_MEMORY);
				command += " " + string(STR_NB_CORES) + " 1";
				command += " " + string(STR_SIMKA_MIN_KMER_SHANNON_INDEX) + " " + _options->getStr(STR_SIMKA_MIN_KMER_SHANNON_INDEX);

				SimkaDistanceParam distanceParams(_options);
				if(distanceParams._computeBrayCurtis) command += " " + STR_SIMKA_DISTANCE_BRAYCURTIS + " ";
				if(distanceParams._computeCanberra) command += " " + STR_SIMKA_DISTANCE_CANBERRA + " ";
				if(distanceParams._computeChord) command += " " + STR_SIMKA_DISTANCE_CHORD + " ";
				if(distanceParams._computeHellinger) command += " " + STR_SIMKA_DISTANCE_HELLINGER + " ";
				if(distanceParams._computeKulczynski) command += " " + STR_SIMKA_DISTANCE_KULCZYNSKI + " ";

	#ifndef CLUSTER
				command += " &";
	#endif

	#ifdef CLUSTER
				string jobFilename = _outputDirTempFilter + "/job_merge/job_merge_" + SimkaAlgorithm<>::toString(i) + ".bash";
				IFile* jobFile = System::file().newFile(jobFilename.c_str(), "w");
				string jobCommand = _jobMergeContents + '\n' + '\n';
				jobCommand += command;

				cout << "\t" << jobCommand << endl;

				jobFile->fwrite(jobCommand.c_str(), jobCommand.size(), 1);
				jobFile->flush();
				string submitCommand = _jobMergeCommand + " " + jobFile->getPath();
				delete jobFile;
				system(submitCommand.c_str());
	#else
				cout << "\t" << command << endl;
				system(command.c_str());
	#endif

				nbJobs += 1;
			}

			if(nbJobs >= _maxJobMerge){
				while(true){
					bool isJobAvailbale = false;

					for(size_t j=0; j<filenameQueue.size(); j++){

						string finishFilename2 = _outputDirTempFilter + "/merge_synchro/" + filenameQueue[j] + ".ok";
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

					if(i >= _bankNames.size()) break;
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
		cout << _nbBanks << endl;


		SimkaDistanceParam distanceParams(_options);
		SimkaStatistics mainStats(_nbBanks, distanceParams);

		for(size_t i=0; i<_nbPartitions; i++){

			Storage* storage = StorageFactory(STORAGE_HDF5).load (_outputDirTempFilter + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".stats");
			LOCAL (storage);

			SimkaStatistics stats(_nbBanks, distanceParams);
			stats.load(storage->getGroup(""));

			cout << stats._nbKmers << endl;
			mainStats += stats;
		}

		mainStats.print();

		mainStats.outputMatrix(_outputDir, _bankNames);
	}


	u_int64_t _maxMemory;
	size_t _nbCores;
	u_int64_t _memoryPerJob;
	size_t _coresPerJob;

	//IBank* _banks;
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
	string _jobCountFilename;
	string _jobMergeFilename;
	string _jobCountCommand;
	string _jobMergeCommand;
	u_int64_t _nbAskedPartitions;

	string _jobCountContents;
	string _jobMergeContents;

	IteratorListener* _progress;
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


class SimkaPotara : public Tool{

public:

	SimkaPotara();
	void execute();
};





#endif
