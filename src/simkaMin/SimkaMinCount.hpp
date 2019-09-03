/*****************************************************************************
 *   SimkaMin: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2019  INRIA
 *   Authors: G.Benoit
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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOUNT_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOUNT_HPP_

/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

//#include "../core/SimkaUtils.hpp"
//#include "Simka2Utils.hpp"
//#include "../minikc/MiniKC.hpp"
//#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
//#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
//#include "../utils/SimkaIoUtils.hpp"
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"

#include "SimkaMinCommons.hpp"
#include "SimkaCommons.hpp"
#include "MurmurHash3.h"
#include <mutex>
//#include "../../thirdparty/KMC/kmc_api/kmc_file.h"
//#include "../../thirdparty/KMC/kmc_api/kmer_defs.h"
//#include "../utils/MurmurHash3.h"

//#define MERGE_BUFFER_SIZE 10000


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
template<size_t span=KMER_DEFAULT_SPAN>
class SelectKmersCommand
{
public:


	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;
    typedef typename Kmer<span>::Type  KmerType;
    typedef typename Kmer<span>::ModelCanonical::Kmer  KmerCanonicalType;


    //typedef typename ModelCanonical::Kmer Lol;

	size_t _kmerSize;
	size_t _sketchSize;
	u_int32_t _seed;
	//vector<u_int64_t> _minHashValues;
	//vector<u_int64_t> _minHashKmers;
	ModelCanonical _model;
	ModelCanonicalIterator _itKmer;

	u_int64_t _hash_otpt[2];

	bool _isMaster;

    //size_t _bufferIndex;
	//size_t _partitionId;

	//vector<u_int64_t> _bufferKmers;
	//vector<u_int32_t> _bufferCounts;

	//vector<u_int64_t> _minHashValues;
	//vector<u_int64_t>& _minHashValuesSynchronized;
	//vector<u_int64_t> _minHashKmers;
	//vector<u_int32_t> _minHashKmersCounts;

	struct KmerCountSorter{
		bool operator() (u_int64_t l, u_int64_t r) { return r > l; }
	};
	//typedef typename KmerCountSorter KmerSorter;

	std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> _kmerCountSorter;

	vector<KmerCanonicalType> _kmers;

	//std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter>& _kmerCountSorterSynch;
	//KmerCountDictionaryType& _kmerCountsSynch;

	Bloom<KmerType>* _bloomFilter;
	u_int64_t _nbInsertedKmersInBloom;

	vector<u_int64_t>& _hashedKmers;
	KmerCountDictionaryType& _kmerCounts;
	//ofstream _outputFile;
	bool _useAbundanceFilter;

	SelectKmersCommand(size_t kmerSize, size_t sketchSize, u_int32_t seed, Bloom<KmerType>* bloomFilter, vector<u_int64_t>& kmers, KmerCountDictionaryType& kmerCounts, bool useAbundanceFilter)
	: _model(kmerSize), _itKmer(_model), _bloomFilter(bloomFilter), _hashedKmers(kmers), _kmerCounts(kmerCounts)
	{
		_kmerSize = kmerSize;
		_sketchSize = sketchSize;
		_seed = seed;
		_isMaster = true;
		_nbInsertedKmersInBloom = 0;
		_useAbundanceFilter = useAbundanceFilter;
	}

	SelectKmersCommand(const SelectKmersCommand& copy)
	: _model(copy._kmerSize), _itKmer(_model), _bloomFilter(copy._bloomFilter), _hashedKmers(copy._hashedKmers), _kmerCounts(copy._kmerCounts)
	{
		_kmerSize = copy._kmerSize;
		_sketchSize = copy._sketchSize;
		_seed = copy._seed;
		_isMaster = false;
		_nbInsertedKmersInBloom = 0;
		_useAbundanceFilter = copy._useAbundanceFilter;
	}


	~SelectKmersCommand(){

		if(_isMaster) return;
		if(_kmerCountSorter.size() == 0) return;
		//cout << "deleteeeeee" << endl;


		size_t sketchSize = _kmerCountSorter.size();
		//cout << sketchSize << endl;
		for(size_t i=0; i<sketchSize; i++){

			//cout << kmers.size()-1-i << endl;
			u_int64_t kmer = _kmerCountSorter.top();
			//cout << kmer << endl;
			//cout << i << ": " << kmer << endl;
			_hashedKmers[_hashedKmers.size()-1-i] = kmer;
			_kmerCountSorter.pop();

		}

	}

	/*
	string revcomp(string& seq){
		string rev = "";
		for(size_t i=seq.size()-1; i>=0; i++){
			if(seq[i] == 'A'){
				rev += 'T';
			}
			else if(seq[i] == 'C'){
				rev += 'G';
			}
			else if(seq[i] == 'G'){
				rev += 'C';
			}
			else if(seq[i] == 'T'){
				rev += 'A';
			}
		}
		return rev;
	}
	*/

	//void minRevComp(string& kmer){
		//string revKmer =
	//}


	void operator()(Sequence& sequence){

		_model.build(sequence.getData(), _kmers);
		//_itKmer.setData(sequence.getData());
		//cout << sequence.toString() << endl;

        //size_t len  = sequence.getDataSize() - _kmerSize + 1; /// _kmerSize;
        //char*  data = sequence.getDataBuffer();

        for(size_t i=0; i<_kmers.size(); i++){

        	KmerCanonicalType& kmer = _kmers[i];
			// We iterate the sequence data by block of size kmerSize
			//for (size_t i=0; i<len; i++, data += 1)
			//{
            // We get the kmer value of the current block
        	//KmerType2 kmer = _model.codeSeed (data, sequence.getDataEncoding());
        	//cout << kmer.value().toString(_kmerSize) << endl;

        	if(!kmer.isValid()) continue;

			//KmerType kmer = _itKmer->value();
			//KmerType kmerRev = revcomp(kmer.value(), _kmerSize);
			//string kmerStr = kmer.value().toString(_kmerSize);
			//string kmerStrRev = kmerRev.toString(_kmerSize);

			//if(kmerStrRev < kmerStr){
			//	kmerStr = kmerStrRev;
			//}

			u_int64_t kmerValue = kmer.value().getVal();
			u_int64_t kmerHashed;
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
			kmerHashed = _hash_otpt[0];

			//cout << kmerStr << ": " << kmerHashed << endl;
			//todo: verifier dabord si le kmer peut etre insérer, plus rapide que els accès au table de hachage (bloom et selected)

			//cout << _useAbundanceFilter << endl;
			if(_useAbundanceFilter){
				processFiltered(kmer.value(), kmerHashed);
			}
			else{
				processUnfiltered(kmerHashed);
			}


            //cout << kmer.isValid() << endl;
            // We update the occurrences number for this kmer value
            //distrib [kmer.value().toInt()] += 1;
        }

    	/*
		for(_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			//cout << _itKmer->value().toString(_kmerSize) << endl;


			Lol kkaka = _itKmer->value().value();

			//cout << _itKmer->value().isValid() << endl;
			KmerType kmer = _itKmer->value();
			KmerType kmerRev = revcomp(kmer, _kmerSize);
			string kmerStr = kmer.toString(_kmerSize);
			string kmerStrRev = kmerRev.toString(_kmerSize);

			if(kmerStrRev < kmerStr){
				kmerStr = kmerStrRev;
			}

			//u_int64_t kmerValue = kmer.getVal();
			u_int64_t kmerHashed;
			MurmurHash3_x64_128 ( kmerStr.c_str(), _kmerSize, 42, &_hash_otpt);
			kmerHashed = _hash_otpt[0];

			if(kmerHashed == 66908235404){
				//cout << kmer.value().isValid() << endl;
				cout << sequence.toString() << endl;
				cout << kmerStr << endl;
			}
			//cout << kmerStr << ": " << kmerHashed << endl;
			//todo: verifier dabord si le kmer peut etre insérer, plus rapide que els accès au table de hachage (bloom et selected)

			//cout << _useAbundanceFilter << endl;
			if(_useAbundanceFilter){
				processFiltered(kmer, kmerHashed);
			}
			else{
				processUnfiltered(kmerHashed);
			}


		}*/
	}

	inline void processUnfiltered(u_int64_t kmerHashed){

		if(_kmerCountSorter.size() < _sketchSize){
			if(_kmerCounts.find(kmerHashed) == _kmerCounts.end()){
				_kmerCountSorter.push(kmerHashed);
				_kmerCounts[kmerHashed] = 1;
				//cout << _kmerCountSorter.size() << endl;
			}
			else{
				_kmerCounts[kmerHashed] += 1;
			}
		}
		else{
			if(kmerHashed < _kmerCountSorter.top()){
				if(_kmerCounts.find(kmerHashed) == _kmerCounts.end()){
					//cout << kmer << "     " << _kmerCounts.size() << endl;
					u_int64_t greaterValue = _kmerCountSorter.top();
					_kmerCounts.erase(greaterValue);
					_kmerCountSorter.pop();
					_kmerCountSorter.push(kmerHashed);
					_kmerCounts[kmerHashed] = 1;
				}
				else{
					_kmerCounts[kmerHashed] += 1;
				}
			}
		}
	}


	inline void processFiltered(const KmerType& kmer, u_int64_t kmerHashed){

		if(_kmerCountSorter.size() < _sketchSize){
			if(_bloomFilter->contains(kmer)){
				//Filling the queue with first elements
				if(_kmerCounts.find(kmerHashed) == _kmerCounts.end()){
					_kmerCountSorter.push(kmerHashed);
					_kmerCounts[kmerHashed] = 2;
					//cout << _kmerCountSorter.size() << endl;
				}
				else{
					_kmerCounts[kmerHashed] += 1;
				}
			}
			else{
				_bloomFilter->insert(kmer);
				_nbInsertedKmersInBloom += 1;
			}
		}
		else{
			if(kmerHashed < _kmerCountSorter.top()){
				if(_bloomFilter->contains(kmer)){

					if(_kmerCounts.find(kmerHashed) == _kmerCounts.end()){
						//cout << kmer << "     " << _kmerCounts.size() << endl;
						u_int64_t greaterValue = _kmerCountSorter.top();
						_kmerCounts.erase(greaterValue);
						_kmerCountSorter.pop();
						_kmerCountSorter.push(kmerHashed);
						_kmerCounts[kmerHashed] = 2;
					}
					else{
						_kmerCounts[kmerHashed] += 1;
					}
				}
				else{
					_bloomFilter->insert(kmer);
					_nbInsertedKmersInBloom += 1;
				}
			}
		}
	}

};
/*
 	size_t _kmerSize;
	size_t _sketchSize;
	vector<u_int64_t> _minHashValues;
	vector<u_int64_t> _minHashKmers;
	ModelCanonical _model;
	ModelCanonicalIterator _itKmer;
	//ModelCanonical model (kmerSize);
	//ModelCanonical::Kmer kmer = model.codeSeed (seq, Data::ASCII);

	pthread_mutex_t* _mutex;
	vector<u_int64_t>& _minHashValuesSynchronized;
	vector<u_int64_t>& _minHashKmersSynchronized;

	MinhashSketcher(size_t kmerSize, size_t sketchSize, pthread_mutex_t* mutex, vector<u_int64_t>& minHashValuesSynchronized, vector<u_int64_t>& minHashKmersSynchronized)
	: _model(kmerSize), _itKmer(_model), _mutex(mutex), _minHashValuesSynchronized(minHashValuesSynchronized), _minHashKmersSynchronized(minHashKmersSynchronized)
	{
		_kmerSize = kmerSize;
		_sketchSize = sketchSize;

		ModelCanonical _model(_kmerSize);
		_minHashValues = vector<u_int64_t>(_sketchSize, -1);
		_minHashKmers = vector<u_int64_t>(_sketchSize, 0);
	}

	MinhashSketcher(const MinhashSketcher& copy)
	: _model(copy._kmerSize), _itKmer(_model), _mutex(copy._mutex), _minHashValuesSynchronized(copy._minHashValuesSynchronized), _minHashKmersSynchronized(copy._minHashKmersSynchronized)
	{
		_kmerSize = copy._kmerSize;
		_sketchSize = copy._sketchSize;
		_minHashValues = vector<u_int64_t>(_sketchSize, -1);
		_minHashKmers = vector<u_int64_t>(_sketchSize, 0);
	}

	~MinhashSketcher(){

		//cout << "deleteeeeee" << endl;
		pthread_mutex_lock(_mutex);

		for(size_t i=0; i<_sketchSize; i++){
			if(_minHashValues[i] < _minHashValuesSynchronized[i]){
				_minHashValuesSynchronized[i] = _minHashValues[i];
				_minHashKmersSynchronized[i] = _minHashKmers[i];
				//cout << _minHashKmers[i] << endl;
			}
		}

		pthread_mutex_unlock(_mutex);
	}



*/




/*
template<size_t span=KMER_DEFAULT_SPAN>
class StorageItKmerCount
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;

    StorageItKmerCount(Iterator<Count>* it){
    	_it = it;
    }

    ~StorageItKmerCount(){
    	delete _it;
    }


	bool next(){
		_it->next();
		return !_it->isDone();
	}

	Count& item(){
		return _it->item();
	}

    Iterator<Count>* _it;
};













template<size_t span=KMER_DEFAULT_SPAN>
class SimkaPartitionWriter
{
public:


    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    //typedef tuple<Count, StorageItKmerCount<span>*> KmerCount_It;
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


	string _outputDir;
	size_t _nbPartitions;

	vector<u_int64_t> _nbKmerPerParts;
	vector<u_int64_t> _nbDistinctKmerPerParts;
	vector<u_int64_t> _chordNiPerParts;
	vector<Bag<Kmer_BankId_Count>* > _bags;
	vector<Bag<Kmer_BankId_Count>* > _cachedBags;


	SimkaPartitionWriter(const string& oututDir, size_t nbPartitions){
		_outputDir = oututDir;
		_nbPartitions = nbPartitions;

		_nbKmerPerParts = vector<u_int64_t>(_nbPartitions, 0);
		_nbDistinctKmerPerParts =  vector<u_int64_t>(_nbPartitions, 0);
		_chordNiPerParts = vector<u_int64_t>(_nbPartitions, 0);


		//vector<Bag<Kmer_BankId_Count>* > bags;
		//vector<Bag<Kmer_BankId_Count>* > cachedBags;
		for(size_t i=0; i<_nbPartitions; i++){
			//string outputFilename = _outputDir + "/" + _datasetID + "_" + Stringify::format("%i", i) + ".gz";
			string outputFilename = _outputDir + "/" + Stringify::format("%i", i) + ".gz";
			Bag<Kmer_BankId_Count>* bag = new BagGzFile<Kmer_BankId_Count>(outputFilename);
			Bag<Kmer_BankId_Count>* cachedBag = new BagCache<Kmer_BankId_Count>(bag, 10000);
			_cachedBags.push_back(cachedBag);
			//BagCache bagCache(*bag, 10000);
			_bags.push_back(bag);
		}
	}

	void insert(u_int64_t kmer, u_int64_t bankId, u_int64_t abundance){


		//kmer.to_long(kmer_bin);
		size_t part = korenXor(kmer) % _nbPartitions; //hash_kmer(kmer_bin) % _nbPartitions;

		Type type; //(kmer_bin[0]);
		//type.setVal(kmer_bin[0]);
		type.setVal(kmer);
		//size_t part = oahash(kmer) % _nbPartitions;
		_cachedBags[part]->insert(Kmer_BankId_Count(type, bankId, abundance));
		_nbDistinctKmerPerParts[part] += 1;
		_nbKmerPerParts[part] += abundance;
		_chordNiPerParts[part] += pow(abundance, 2);

	}

	void end(){
		for(size_t i=0; i<_nbPartitions; i++){
			//bags[i]->flush();
			//cachedBags[i]->flush();
			delete _cachedBags[i];
			//delete bags[i];
		}


		for(size_t i=0; i<_nbPartitions; i++){
			string outputFilename = _outputDir + "/" + Stringify::format("%i", i) + ".gz";
			checkGzFile(outputFilename);
		}
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


	inline u_int64_t korenXor(u_int64_t x)const{
	        x ^= (x << 21);
	        x ^= (x >> 35);
	        x ^= (x << 4);
	        return x;
	}
};
*/




/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2ComputeKmerSpectrumAlgorithm : public Algorithm
{

public:

	typedef typename Kmer<span>::Type  KmerType;
	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;
	struct KmerCountSorter{
		bool operator() (u_int64_t l, u_int64_t r) { return r > l; }
	};

	//struct kxpcomp { bool operator() (KmerCount_It& l,KmerCount_It& r) { return (r._count.value < l._count.value); } } ;

	//u_int64_t _nbReads;

	//size_t _nbPartitions;
	u_int64_t _maxMemory;
	size_t _nbCores;
	string _outputDir;
	string _outputDirTemp;
	//size_t _nbBanks;
	string _inputFilename;
	//string _datasetID;
	u_int8_t _kmerSize;
	//pair<CountNumber, CountNumber> _abundanceThreshold;
	//SIMKA_SOLID_KIND _solidKind;
	//bool _soliditySingle;
	int64_t _maxNbReads;
	size_t _minReadSize;
	double _minReadShannonIndex;
	//double _minKmerShannonIndex;
	//size_t _nbMinimizers;
	//size_t _nbCores;

	//SimkaStatistics* _stats;
	//SimkaDistance* _simkaDistance;

	//string _banksInputFilename;
	//string _h5Filename;
	//vector<string> _tempFilenamesToDelete;
	//IBank* _banks;
	IProperties* _options;
	//size_t _localNbPartitions;

	//vector<string> _bankNames;
	//vector<u_int64_t> _nbReadsPerDataset;

	//string _outputFilenameSuffix;

	//u_int64_t _totalKmers;
	//vector<size_t> _nbBankPerDataset;
	//size_t _nbBankPerDataset;

	//string _largerBankId;
	//bool _computeSimpleDistances;
	//bool _computeComplexDistances;
	//bool _keepTmpFiles;

	//string _kmerDatataseFilename;

	//vector<ICommand*> _cmds;
	//SimkaPartitionWriter<span>* _partitionWriter;



	u_int32_t _seed;
	u_int32_t _sketchSize;
	bool _useAbundanceFilter;
	u_int32_t _nbDatasets;
	//pthread_mutex_t _mutex;

	//typedef typename SelectKmersCommand<span>::KmerCountSorter KmerCountSorter;
	//std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> _kmerCountSorter;
	//KmerCountDictionaryType _kmerCounts;

	//size_t _nbBanks;
	//vector<string> _bankNames;
	//vector<size_t> _nbBankPerDataset;

	vector<thread*> _threads;
	size_t _maxRunningThreads;
	vector<size_t> _runningThreadIds;
	size_t _nbRunningThreads;
	vector<size_t> _finishedThreads;
	mutex countKmersMutex;

	ofstream _outputFile;
	//string _outputFilenameKmers;
	//string _outputFilenameIds;


	IteratorListener* _progress;
	u_int64_t _progress_nbDatasetsToProcess;
	u_int64_t _progress_nbDatasetsProcessed;
	string _progress_text;

	Simka2ComputeKmerSpectrumAlgorithm(IProperties* options):
		Algorithm("simka", -1, options)
	{
	}

	void execute(){
		//pthread_mutex_init(&_mutex, NULL);


		parseArgs();
		createDirs();

		cout << endl << "Checking input file validity..." << endl;
		SimkaCommons::checkInputValidity(_outputDirTemp, _inputFilename, _progress_nbDatasetsToProcess);

		_progress = this->createIteratorListener (_progress_nbDatasetsToProcess, ""); //new ProgressSynchro (
			//this->createIteratorListener (_progress_nbDatasetsToProcess, ""),
			//System::thread().newSynchronizer());
	    _progress->setMessage (Stringify::format (_progress_text.c_str(), _progress_nbDatasetsProcessed, _progress_nbDatasetsToProcess));
		_progress->init ();

		countDatasets();


		string command = "rm -rf " + _outputDirTemp;
		system(command.c_str());

		cout << "Output results: " << _outputDir << endl;
	}

	void parseArgs(){

		_options = getInput();

		_seed = _options->getInt(STR_SIMKA_SEED);
		_sketchSize = _options->getInt(STR_SIMKA_SKETCH_SIZE);
		_useAbundanceFilter = _options->get(STR_SIMKA_ABUNDANCE_FILTER);

		_maxMemory = _options->getInt(STR_MAX_MEMORY);
	    _nbCores = _options->getInt(STR_NB_CORES);
		_inputFilename = _options->getStr(STR_URI_INPUT);
		//_datasetID = _options->getStr(STR_SIMKA2_DATASET_ID);
		_outputDir = _options->getStr(STR_URI_OUTPUT); // ? _options->getStr(STR_URI_OUTPUT) : "./";
		if(_outputDir.empty()) _outputDir = "./simkaMin_kmers.bin";
		_outputDirTemp = System::file().getDirectory(_outputDir) + "/__simkaMin_temp__/";
		//cout << "outputdir temp to check: " << _outputDirTemp << endl;
		//_outputDirTemp = _options->get(STR_URI_OUTPUT_TMP) ? _options->getStr(STR_URI_OUTPUT_TMP) : "./";
		_kmerSize = _options->getInt(STR_KMER_SIZE);
		//_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
		//_abundanceThreshold.second = min((u_int64_t)_options->getInt(STR_KMER_ABUNDANCE_MAX), (u_int64_t)(999999999));

		//_nbPartitions = _options->getInt(STR_SIMKA2_NB_PARTITION);
		//cout << _options->getInt(STR_KMER_ABUNDANCE_MAX) << endl;
		//cout << _abundanceThreshold.second << endl;
		//_soliditySingle = _options->get(STR_SIMKA_SOLIDITY_PER_DATASET);
		//_nbMinimizers = _options->getInt(STR_KMER_PER_READ);
		//_maxDisk = getInput()->getInt(STR_MAX_DISK);

		//read filter
		_maxNbReads = _options->getInt(STR_SIMKA_MAX_READS);
		_minReadSize = _options->getInt(STR_SIMKA_MIN_READ_SIZE);
		_minReadShannonIndex = _options->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
		_minReadShannonIndex = std::max(_minReadShannonIndex, 0.0);
		_minReadShannonIndex = std::min(_minReadShannonIndex, 2.0);


		if(!System::file().doesExist(_inputFilename)){
			std::cerr << "Error: input does not exist (" << _inputFilename << ")" << std::endl;
			exit(1);
		}

		if(System::file().doesExist(_outputDir)){
			std::cerr << "Error: output file already exist (" << _outputDir << ")" << std::endl;
			exit(1);
		}


		_progress_text = "Sketching datasets (%d/%d)";
		//_nbBankPerDataset = _options->getInt("-nb-dataset");

		//_minKmerShannonIndex = _options->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
		//_minKmerShannonIndex = std::max(_minKmerShannonIndex, 0.0);
		//_minKmerShannonIndex = std::min(_minKmerShannonIndex, 2.0);

		//if(!System::file().doesExist(_inputFilename)){
		//	cerr << "ERROR: Input filename does not exist" << endl;
		//	exit(1);
		//}

		//if(!System::file().doesExist(_outputDir)){
		//	std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
		//	exit(EXIT_FAILURE);
			/*
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
				std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
				exit(1);
			}*/
		//}

		//_outputDirTemp = _outputDirTemp;

		//if(!System::file().doesExist(_outputDirTemp)){
		//std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
		//exit(EXIT_FAILURE);
			/*
			int ok = System::file().mkdir(_outputDirTemp, -1);
			if(ok != 0){
				std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
				exit(1);
			}*/
		//}

		//_outputDirTemp = System::file().getRealPath(_outputDirTemp) + "/";
		//cout << _outputDirTemp << endl;
		//_outputDirTemp += "/" + _datasetID + "_temp" + "/";
		//System::file().mkdir(_outputDirTemp, -1);

		//_options->setStr(STR_URI_OUTPUT_TMP, _outputDirTemp);
		//System::file().mkdir(_outputDirTemp + "/input/", -1);

		//_maxMemory = _maxMemory / 1000;
		//_maxMemory = max(_maxMemory, (u_int64_t) 1);

		/*
		if(_outputDir.empty()){
			_outputDir = "./simkaMin_kmers.bin";
		}
		else if (_outputDir.find(".") == std::string::npos){
			_outputDir += ".bin";
		}
		_outputDir = System::file().getBaseName(_outputDir);

		cout << endl << endl;
		cout << _outputDir << endl;
		vector<string> fields;
		stringstream outputFilenameStream(_outputDir);
		string field;

		while(std::getline(outputFilenameStream, field, '.'))
		{
			cout << field << endl;
			fields.push_back(field);
		}

		string prefix = fields[0];
		string extension = "";
		for(size_t i=1; i<fields.size(); i++){
			extension += fields[i];
		}*/

		//_outputFilenameIds = _outputDir + ".ids";
		//_outputFilenameKmers = _outputDir + ".kmers";

		//System::file().remove(_outputFilenameIds);
		//System::file().remove(_outputFilenameKmers);

	}


	void createDirs(){

		//if(!System::file().doesExist(_outputDir)){
		//int ok = System::file().mkdir(_outputDir, -1);
		//if(ok != 0){
		//      std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
		//      exit(1);
		//}
		//}

		if(!System::file().doesExist(_outputDirTemp)){
			int ok = System::file().mkdir(_outputDirTemp, -1);
			if(ok != 0){
		        std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
		        exit(1);
			}
		}

		_outputDirTemp = System::file().getRealPath(_outputDirTemp) + "/";



		//_outputDirTemp += "/simka_output_temp/";
		//System::file().mkdir(_outputDirTemp, -1);

		//_args->setStr(STR_URI_OUTPUT_TMP, _outputDirTemp);
		//System::file().mkdir(_outputDirTemp + "/input/", -1);

	}


	void countDatasets(){

		//cout << endl << endl;
		//cout << "Sketching..." << endl;


		_outputFile.open(_outputDir, ios::binary);

		//Save sketch info
		//u_int8_t kmerSize = _kmerSize;
		//u_int32_t sketchSize = _sketchSize;
		//u_int32_t seed = _seed;
		_nbDatasets = 0;
		_outputFile.write((const char*)&_kmerSize, sizeof(_kmerSize));
		_outputFile.write((const char*)&_sketchSize, sizeof(_sketchSize));
		_outputFile.write((const char*)&_seed, sizeof(_seed));
		_outputFile.write((const char*)&_nbDatasets, sizeof(_nbDatasets));

		//cout << _maxRunningThreads << endl;

		size_t threadId = 0;
		//vector<thread> threads; //(_nbCores);
		//_isThreadRunning = vector<bool>(_nbCores);
		_nbRunningThreads = 0;
		_maxRunningThreads = _nbCores;









		string inputDir = _outputDirTemp; // + "/input/";
		ifstream inputFile(_inputFilename.c_str());

		//ofstream outputFileIds(_outputFilenameIds.c_str(), ios::binary);
		//_banksInputFilename =  inputDir + "__input_simka__"; //_inputFilename + "_dsk_dataset_temp__";
		//IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		string line;
		string linePart;
		vector<string> lineIdDatasets;
		vector<string> linepartPairedDatasets;
		vector<string> linepartDatasets;

		//string bankFileContents = "";

		size_t datasetId = 0;
		u_int64_t lineIndex = 0;
		u_int64_t bankIdBytePos = 0;

		while(getline(inputFile, line)){

			line.erase(std::remove(line.begin(),line.end(),' '),line.end());
			if(line == "") continue;

			//cout << line << endl;
			lineIdDatasets.clear();
			linepartPairedDatasets.clear();
			//vector<string> filenames;

			stringstream lineStream(line);
			while(getline(lineStream, linePart, ':')){
				lineIdDatasets.push_back(linePart);
			}

			string bankId = lineIdDatasets[0];
			string linePairedDatasets = lineIdDatasets[1];

			stringstream linePairedDatasetsStream(linePairedDatasets);
			while(getline(linePairedDatasetsStream, linePart, ';')){
				linepartPairedDatasets.push_back(linePart);
			}

			string subBankFilename = inputDir + bankId;
			IFile* subBankFile = System::file().newFile(subBankFilename, "wb");
			//cout << subBankFile->getPath() << endl;
			string subBankContents = "";
			size_t nbBankPerDataset = linepartPairedDatasets.size();

			for(size_t i=0; i<linepartPairedDatasets.size(); i++){
				string lineDatasets = linepartPairedDatasets[i];

				linepartDatasets.clear();

				stringstream lineDatasetsStream(lineDatasets);
				while(getline(lineDatasetsStream, linePart, ',')){
					linepartDatasets.push_back(linePart);
					//cout << "\t" << linePart << endl;
				}

				//bankFileContents += linepartDatasets[0] + "\n";


				for(size_t i=0; i<linepartDatasets.size(); i++){
					string filename = linepartDatasets[i];
					if(filename.at(0) == '/'){
						subBankContents +=  filename + "\n";
					}
					else{
						string dir = System::file().getRealPath(_inputFilename);
						dir = System::file().getDirectory(dir);
						subBankContents +=  dir + "/" + filename + "\n";
					}
				}

			}

			subBankContents.erase(subBankContents.size()-1);
			subBankFile->fwrite(subBankContents.c_str(), subBankContents.size(), 1);
			subBankFile->flush();
			delete subBankFile;

			//bankFileContents += inputDir + "/" + bankId + "\n";
			lineIndex += 1;

			startNewThread(datasetId, subBankFilename, nbBankPerDataset);
			//count();
			//_bankNames.push_back(bankId);

			datasetId += 1;
			_nbDatasets += 1;
		}



		//bankFileContents.erase(bankFileContents.size()-1);
		//bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);
		//bankFile->flush();
		//delete bankFile;


		joinThreads();
		_progress->finish();

		inputFile.close();

		writeIds();

		//outputFileIds.seekp(0);
		//outputFileIds.write((const char*)&nbDatasets, sizeof(nbDatasets));
		_outputFile.close();
		//outputFileIds.close();

	}

	void writeIds(){

		_outputFile.seekp(SimkaMinCommons::getFilePosition_nbDatasets());
		_outputFile.write((const char*)&_nbDatasets, sizeof(_nbDatasets));
		_outputFile.seekp(SimkaMinCommons::getFilePosition_sketchIds(_nbDatasets, _sketchSize));

		ifstream inputFile(_inputFilename.c_str());

		string line;
		string linePart;
		vector<string> lineIdDatasets;


		while(getline(inputFile, line)){

			line.erase(std::remove(line.begin(),line.end(),' '),line.end());
			if(line == "") continue;

			lineIdDatasets.clear();

			stringstream lineStream(line);
			while(getline(lineStream, linePart, ':')){
				lineIdDatasets.push_back(linePart);
			}

			string bankId = lineIdDatasets[0];

			u_int8_t idSize = bankId.size();
			_outputFile.write((const char*)& idSize, sizeof(idSize));
			_outputFile.write(bankId.c_str(), bankId.size());

		}

		inputFile.close();

	}


	void startNewThread(size_t datasetId, const string& inputFilename, size_t nbBankPerDataset){



		//for (size_t i=0; i<_nbBanks; i++){
		//	cout << i << endl;
		thread* t = new thread(&Simka2ComputeKmerSpectrumAlgorithm::countKmersOfDataset, this, datasetId, inputFilename, nbBankPerDataset);
		_threads.push_back(t);
		_runningThreadIds.push_back(datasetId);
		//threadId += 1;
		_nbRunningThreads += 1;
		//_isThreadRunning[threadId] = true;
		//_nbRunningThreads[i] += 1;

		if(_nbRunningThreads >= _maxRunningThreads){
			waitThreads();
		}


		//}

		//string filename = _outputDirTemp + "/selectedKmers.bin";
		//ofstream selectKmersFile(filename.c_str(), ios::binary);

		//cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;
		//_selectedKmerSorter.pop(); //there is always one extra element because of a >= optimization...
		//cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;
		//u_int64_t size = _selectedKmerSorter.size();
		//for(size_t i=0; i<size; i++){
		//	u_int64_t kmerValue = _selectedKmerSorter.top();
		//	_selectedKmerSorter.pop();
		//	_selectedKmersIndex[kmerValue] = i;
		//	//selectKmersFile.write((const char*)&kmerValue, sizeof(kmerValue)); //todo mphf loading can be done in memory, not required to write all selected kmers on disk
		//}
		//_nbUsedKmers = _selectedKmersIndex.size() ;
		//cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;

		//selectKmersFile.close();
	}

	//vector<unordered_map<u_int64_t, KmerCountType> > _skecthCounts;

	//unordered_map<u_int64_t, vector<KmerCountType> > _;

	void countKmersOfDataset(size_t datasetId, const string& inputFilename, size_t nbBankPerDataset){

		//TODO lock probably not required
		//countKmersMutex.lock();
		//cout << "start: " << inputFilename << endl;
		//countKmersMutex.unlock();

		IBank* bank = Bank::open(inputFilename);
		LOCAL(bank);

		SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);
		IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _maxNbReads, nbBankPerDataset);

		LOCAL(filteredBank);

		Iterator<Sequence>* itSeq = filteredBank->iterator();
		LOCAL(itSeq);
		//Iterator<Sequence>* itSeq = createIterator<Sequence> (
		//		filteredBank->iterator(),
		//		filteredBank->estimateNbItems(),
		//		"Computing minhash sketch and counting"
		//);
		//LOCAL(itSeq);


		IDispatcher* dispatcher = new SerialDispatcher();
		Bloom<KmerType>*  bloomFilter = 0;

		if(_useAbundanceFilter){
			u_int64_t bloomMemoryBits = (_maxMemory * MBYTE * 8)  / _maxRunningThreads;
			bloomMemoryBits = max(bloomMemoryBits, (u_int64_t) 10000);
			bloomFilter = new BloomCacheCoherent<KmerType>(bloomMemoryBits, 7);
		}
		//mutex commandMutex;
		//std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> kmerCountSorter;
		//unordered_map<u_int64_t, KmerCountType> kmerCounts;


		vector<u_int64_t> kmers(_sketchSize, 0); //TODO only used for reversing kmers not really optimized...
		KmerCountDictionaryType _kmerCounts;

		{
			SelectKmersCommand<span> command(_kmerSize, _sketchSize, _seed, bloomFilter, kmers, _kmerCounts, _useAbundanceFilter);
			dispatcher->iterate (itSeq, command, 1000);
		}

		/*
		ModelCanonical model;
		ModelCanonicalIterator itKmer(model);

		u_int64_t _hash_otpt[2];



		u_int64_t _nbInsertedKmersInBloom = 0;



		for(itSeq->first(); !itSeq->isDone(); itSeq->next()){

			Sequence& sequence = itSeq->item();


		}
		*/





		delete dispatcher;
		delete bloomFilter;



		countKmersMutex.lock();















		u_int64_t filePos = (datasetId * _sketchSize * sizeof(KmerAndCountType)) + KMER_SPECTRUM_HEADER_SIZE;
		//cout << "DATASTE ID: " << datasetId << "    " << filePos << endl;
		_outputFile.seekp(filePos);


		//_kmerCountSorter.pop(); //Discard greater element because queue size is always equal to (_sketchSize + 1) because of an optimization

		//cout << "----------" << endl;
		for(size_t i=0; i<kmers.size(); i++){

			u_int64_t kmer = kmers[i];

			//cout << kmer << endl;
			KmerAndCountType kmerCount(kmer, _kmerCounts[kmer]);
			//cout << kmer << endl;
			//KmerCountType count = ;

			_outputFile.write((const char*)&kmerCount, sizeof(kmerCount));
			//_outputFile.write((const char*)&count, sizeof(count));

			//mergeSynch(kmer, _kmerCounts[kmer]);
			//KmerCount kmerCount = _kmerCountSorter.top();
			//cout << kmerCount._kmer << "  " << kmerCount._count << endl;
			//_partitionWriter->insert(kmerCount._kmer, _datasetIDbin, kmerCount._count);
			//_kmerCountSorter.pop();
			//_partitionWriter->insert(_minHashKmers[i], _datasetIDbin, _minHashKmersCounts[i] );
			//cout << _minHashKmers[i] << " " << _minHashKmersCounts[i] << endl;
		}

		System::file().remove(inputFilename);


		_progress_nbDatasetsProcessed += 1;
	    _progress->setMessage (Stringify::format (_progress_text.c_str(), _progress_nbDatasetsProcessed, _progress_nbDatasetsToProcess));
	    _progress->inc(1);

		//cout << "end: " << inputFilename << endl;
		_finishedThreads.push_back(datasetId);

		countKmersMutex.unlock();
	}







	void waitThreads(){
		while(1){



			bool isThreadAvailbale = false;

			countKmersMutex.lock();

			for(size_t i=0; i<_finishedThreads.size(); i++){
				size_t threadId = _finishedThreads[i];

				//_runningThreadIds.erase(std::remove(_runningThreadIds.begin(), _runningThreadIds.end(), threadId), _runningThreadIds.end());
				auto it = find(_runningThreadIds.begin(), _runningThreadIds.end(), threadId);
				int pos = distance(_runningThreadIds.begin(), it);

				//cout << "\t removing thread " <<  threadId << " (pos: "  << pos << ")" << endl;

				_runningThreadIds.erase(_runningThreadIds.begin()+pos);
				_threads[pos]->join();
				delete _threads[pos];
				_threads.erase(_threads.begin()+pos);

				_nbRunningThreads -= 1;
				isThreadAvailbale = true;

			}

			_finishedThreads.clear();

			countKmersMutex.unlock();

			if(isThreadAvailbale){
				//cout << _runningThreadIds.size() << " " << _threads.size() << endl;
				//countKmersMutex.unlock();
				break;
			}


			sleep(1);

		}

	}

	void joinThreads(){
		while(_nbRunningThreads > 0)
			waitThreads();
	}






};











class Simka2ComputeKmerSpectrum : public Tool{
public:


	Simka2ComputeKmerSpectrum(): Tool ("SimkaMin-ComputeKmerSpectrum"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

		//Main parser
	    //parser->push_front (new OptionNoParam (STR_SIMKA_COMPUTE_DATA_INFO, "compute (and display) information before running Simka, such as the number of reads per dataset", false));
	    //parser->push_front (new OptionNoParam (STR_SIMKA_KEEP_TMP_FILES, "keep temporary files", false));
	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_SEED, "seed used for random k-mer selection", false, "100"));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output filename for kmer spectrum", false, "./simkaMin_kmers.bin"));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input filename | TODO SPECIF", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_DATASET_ID, "identifier of the input dataset", true));


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
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "21"));
	    kmerParser->push_back (new OptionOneParam (STR_SIMKA_SKETCH_SIZE, "number of kmers used to compute distances", false, "100000"));
	    kmerParser->push_back (new OptionNoParam (STR_SIMKA_ABUNDANCE_FILTER, "filter out k-mer seen one time (potentially erroneous)", false));
	    //kmerParser->push_back(dskParser->getParser (STR_KMER_SIZE));
	    //kmerParser->push_back(new OptionOneParam (STR_KMER_PER_READ.c_str(), "number of selected kmers per read", false, "0"));
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "1"));
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "2"));
	    //KmerCountType maxAbundance = -1;
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "max abundance a kmer can have to be considered", false, Stringify::format("%i", maxAbundance)));

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
	    //readParser->push_back (new OptionOneParam ("-nb-dataset", "nb paired datasets", true));

	    //Core parser
	    IOptionsParser* coreParser = new OptionsParser ("core");
	    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
	    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (MB). Only used if -filter is enabled", false, "8000"));
	    //coreParser->push_back (new OptionOneParam (STR_SIMKA2_NB_PARTITION, "nb partitions", true));


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


		//IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();

	    //if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE)))  {  p->setDefaultValue ("7"); }
		//parser->push_back(dskParser);
		//if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE))) { p->setDefaultValue ("7"); }
		//dskParser->setVisible(false);

		parser->getParser(STR_NB_CORES)->setVisible(false);
		//getParser()->push_back(parser);
	    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

	    //return parser;
	}

	~Simka2ComputeKmerSpectrum(){

	}

	struct Parameter
	{
	    //Parameter (Simka& simka, IProperties* props) : props(props) {}
	    Parameter (IProperties* props) : _props(props) {}
	    //Simka& _simka;
	    IProperties* _props;
	};

	template<size_t span> struct Functor {

    	void operator ()  (Parameter p){

    		Simka2ComputeKmerSpectrumAlgorithm<span>* algo = new Simka2ComputeKmerSpectrumAlgorithm<span>(p._props);
    		algo->execute();
    		delete algo;
    	}


	};

	void execute ()
	{
		IProperties* input = getInput();
		//Parameter params(*this, getInput());
		Parameter params(input);

		size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

	    Integer::apply<Functor,Parameter> (kmerSize, params);
	}

};








#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOUNT_HPP_ */
