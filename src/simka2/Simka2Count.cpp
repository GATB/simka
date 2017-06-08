/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
#include "../core/SimkaUtils.hpp"
#include "Simka2Utils.hpp"
//#include "../minikc/MiniKC.hpp"
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include "../utils/SimkaIoUtils.hpp"
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"

//#include "../../thirdparty/KMC/kmc_api/kmc_file.h"
//#include "../../thirdparty/KMC/kmc_api/kmer_defs.h"
#include <unordered_map>
#include "../utils/MurmurHash3.h"

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


    typedef typename Kmer<span>::Type  KmerType;
	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;

	size_t _kmerSize;
	size_t _sketchSize;
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
	KmerCountDictionaryType _kmerCounts;

	std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter>& _kmerCountSorterSynch;
	KmerCountDictionaryType& _kmerCountsSynch;

	Bloom<KmerType>* _bloomFilter;
	u_int64_t _nbInsertedKmersInBloom;

	SelectKmersCommand(size_t kmerSize, size_t sketchSize, std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter>& kmerCountSorterSynch, KmerCountDictionaryType& kmerCountsSynch, Bloom<KmerType>* bloomFilter)
	: _model(kmerSize), _itKmer(_model), _kmerCountSorterSynch(kmerCountSorterSynch), _kmerCountsSynch(kmerCountsSynch), _bloomFilter(bloomFilter)
	{
		_kmerSize = kmerSize;
		_sketchSize = sketchSize;
		_isMaster = true;
		_nbInsertedKmersInBloom = 0;
	}

	SelectKmersCommand(const SelectKmersCommand& copy)
	: _model(copy._kmerSize), _itKmer(_model), _kmerCountSorterSynch(copy._kmerCountSorterSynch), _kmerCountsSynch(copy._kmerCountsSynch), _bloomFilter(copy._bloomFilter)
	{
		_kmerSize = copy._kmerSize;
		_sketchSize = copy._sketchSize;
		_isMaster = false;
		_nbInsertedKmersInBloom = 0;
	}


	~SelectKmersCommand(){

		if(_isMaster) return;
		if(_kmerCountSorter.size() == 0) return;
		//cout << "deleteeeeee" << endl;
		/*
		//cout << "deleteeeeee" << endl;
		pthread_mutex_lock(_mutex);

		for(size_t i=0; i<_sketchSize; i++){
			if(_minHashValues[i] < _minHashValuesSynchronized[i]){
				_minHashValuesSynchronized[i] = _minHashValues[i];
				_minHashKmersSynchronized[i] = _minHashKmers[i];
				//cout << _minHashKmers[i] << endl;
			}
		}

		pthread_mutex_unlock(_mutex);*/



		//_mutex.lock();

		//cout << "Inserted k-mers: " << _nbInsertedKmersInBloom << endl;
		//cout << "deleteeeeee" << endl;
		//cout << this << endl;

		_kmerCountSorter.pop(); //Discard greater element because queue size is always equal to (_sketchSize + 1) because of an optimization

		size_t sketchSize = _kmerCountSorter.size();
		for(size_t i=0; i<sketchSize; i++){
			u_int64_t kmer = _kmerCountSorter.top();
			_kmerCountSorter.pop();
			mergeSynch(kmer, _kmerCounts[kmer]);
			//KmerCount kmerCount = _kmerCountSorter.top();
			//cout << kmerCount._kmer << "  " << kmerCount._count << endl;
			//_partitionWriter->insert(kmerCount._kmer, _datasetIDbin, kmerCount._count);
			//_kmerCountSorter.pop();
			//_partitionWriter->insert(_minHashKmers[i], _datasetIDbin, _minHashKmersCounts[i] );
			//cout << _minHashKmers[i] << " " << _minHashKmersCounts[i] << endl;
		}

		//cout << "deleteeeeee1" << endl;

		//_mutex.unlock();
	}

	inline void mergeSynch(u_int64_t kmer, KmerCountType count){
		if(_kmerCountsSynch.find(kmer) == _kmerCountsSynch.end()){
			if(_kmerCountSorterSynch.size() > _sketchSize){
				if(kmer < _kmerCountSorterSynch.top() ){
					//cout << kmer << "     " << _kmerCounts.size() << endl;
					u_int64_t greaterValue = _kmerCountSorterSynch.top();
					_kmerCountsSynch.erase(greaterValue);
					_kmerCountSorterSynch.pop();
					_kmerCountSorterSynch.push(kmer);
					_kmerCountsSynch[kmer] = count;
				}
				//else{
				//	cout << "\t\tnonon" << endl;
				//}
			}
			else{ //Filling the queue with first elements
				_kmerCountSorterSynch.push(kmer);
				_kmerCountsSynch[kmer] = count;
			}
		}
		else{
			_kmerCountsSynch[kmer] += count;
		}
	}

	void operator()(Sequence& sequence){
		//cout << sequence.getIndex() << endl;
		//cout << sequence.toString() << endl;

		//cout << _isMaster << endl;
		_itKmer.setData(sequence.getData());

		for(_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			//cout << _itKmer->value().toString(_kmerSize) << endl;

			KmerType kmer = _itKmer->value();



			u_int64_t kmerValue = kmer.getVal();
			u_int64_t kmerHashed;
			MurmurHash3_x64_128 ( (const char*)&kmerValue, sizeof(kmerValue), 100, &_hash_otpt);
			kmerHashed = _hash_otpt[0];

			//todo: verifier dabord si le kmer peut etre insérer, plus rapide que els accès au table de hachage (bloom et selected)

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
						if(kmerHashed < _kmerCountSorter.top() ){
							_bloomFilter->insert(kmer);
							_nbInsertedKmersInBloom += 1;
						}
					}
				}


				//else{
					//cout << "test" << endl;
					//}
			}
			//else{
				//if(_kmerCountSorter.size() < _sketchSize){
				//	_bloomFilter->insert(kmer);
				//	_nbInsertedKmersInBloom += 1;
					//cout << "filling1" << endl;
				//}
				//else{
				//	if(kmerHashed < _kmerCountSorter.top() ){
				//		_bloomFilter->insert(kmer);
				//		_nbInsertedKmersInBloom += 1;
				//	}
					//}
			//}


			//u_int64_t kmer = hash_otpt[0];
			//cout << kmer << endl;


			//if(kmer < 100) continue;


	    	//else{
	    	//	_kmerCounts[kmer] += 1;
	    	//}

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

template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:


	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads, size_t nbDatasets) : BankDelegate (ref), _ref2(0), _filter(filter)  {
		_maxReads = maxReads;
		_nbDatasets = nbDatasets;
		setRef2(_ref->iterator ());
	}



	~SimkaPotaraBankFiltered(){

	    std::vector<Iterator<Sequence>*> itBanks =  _ref2->getComposition();
	    for(size_t i=0; i<itBanks.size(); i++){
	    	delete itBanks[i];
	    }

	    //_ref2->
		setRef2(0);
	}

    Iterator<Sequence>* iterator ()
    {
        return new SimkaInputIterator<Sequence, Filter> (_ref2, _nbDatasets, _maxReads, _filter);

    }

private:


    Iterator<Sequence>* _ref2;
    void setRef2 (Iterator<Sequence>* ref2)  { SP_SETATTR(ref2); }

    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadToProcess;
    size_t _datasetId;
    size_t _nbDatasets;
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

	/*
	inline u_int64_t hash_kmer(const vector<uint64>& kmer_bin){
		uint64 result = 0;

	    //LargeInt<precision> intermediate = elem;
	    for (size_t i=0;i<kmer_bin.size();i++)
	    {
	        //chunk = (intermediate & mask).value[0];
	        //intermediate = intermediate >> 64;
	        result ^= korenXor(kmer_bin[i]);
	    }
	    return result;
	}*/

	inline u_int64_t korenXor(u_int64_t x)const{
	        x ^= (x << 21);
	        x ^= (x >> 35);
	        x ^= (x << 4);
	        return x;
	}
};





/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2ComputeKmerSpectrumAlgorithm : public Algorithm
{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    //typedef tuple<Count, StorageItKmerCount<span>*> KmerCount_It;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    struct KmerCount_It{
    	Count _count;
    	StorageItKmerCount<span>* _it;

    	KmerCount_It(){

    	}

    	KmerCount_It(Count& count, StorageItKmerCount<span>* it){
    		_count = count;
    		_it = it;
    	}
    };

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

	//struct kxpcomp { bool operator() (KmerCount_It& l,KmerCount_It& r) { return (r._count.value < l._count.value); } } ;

	u_int64_t _nbReads;

	string _binDir;
	size_t _nbPartitions;
	u_int64_t _maxMemory;
	size_t _nbCores;
	string _outputDir;
	string _outputDirTemp;
	//size_t _nbBanks;
	string _inputFilename;
	string _datasetID;
	size_t _kmerSize;
	pair<CountNumber, CountNumber> _abundanceThreshold;
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
	string _h5Filename;
	//vector<string> _tempFilenamesToDelete;
	//IBank* _banks;
	IProperties* _options;
	size_t _localNbPartitions;

	u_int64_t _datasetIDbin;
	//vector<string> _bankNames;
	//vector<u_int64_t> _nbReadsPerDataset;

	//string _outputFilenameSuffix;

	//u_int64_t _totalKmers;
	//vector<size_t> _nbBankPerDataset;
	size_t _nbBankPerDataset;

	string _largerBankId;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	bool _keepTmpFiles;

	string _kmerDatataseFilename;

	vector<ICommand*> _cmds;
	SimkaPartitionWriter<span>* _partitionWriter;


	vector<vector<u_int64_t>> _bufferKmers;
	vector<vector<u_int32_t>> _bufferCounts;
	vector<size_t> _bufferIndex;

	vector<u_int64_t> _minHashValues;
	vector<u_int64_t> _minHashKmers;
	vector<u_int32_t> _minHashKmersCounts;

	size_t _sketchSize;
	pthread_mutex_t _mutex;

	typedef typename SelectKmersCommand<span>::KmerCountSorter KmerCountSorter;
	std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> _kmerCountSorter;
	KmerCountDictionaryType _kmerCounts;


	Simka2ComputeKmerSpectrumAlgorithm(IProperties* options, const string& execFilename):
		Algorithm("simka", -1, options)
	{
		_binDir = System::file().getDirectory(execFilename);
	}

	void execute(){
		_datasetIDbin = 0;
		pthread_mutex_init(&_mutex, NULL);

		parseArgs();
		//layoutInputFilename();
		count();
		partitionKmerCounts();
		saveInfos();


	}

	void parseArgs(){

		_options = getInput();

		_sketchSize = _options->getInt(STR_SIMKA_SKETCH_SIZE);

		_computeSimpleDistances = _options->get(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
		_computeComplexDistances = _options->get(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);
		_keepTmpFiles = _options->get(STR_SIMKA_KEEP_TMP_FILES);
		_maxMemory = _options->getInt(STR_MAX_MEMORY);
	    _nbCores = _options->getInt(STR_NB_CORES);
		_inputFilename = _options->getStr(STR_URI_INPUT);
		_datasetID = _options->getStr(STR_SIMKA2_DATASET_ID);
		_outputDir = _options->get(STR_URI_OUTPUT) ? _options->getStr(STR_URI_OUTPUT) : "./";
		_outputDirTemp = _options->get(STR_URI_OUTPUT_TMP) ? _options->getStr(STR_URI_OUTPUT_TMP) : "./";
		_kmerSize = _options->getInt(STR_KMER_SIZE);
		_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
		_abundanceThreshold.second = min((u_int64_t)_options->getInt(STR_KMER_ABUNDANCE_MAX), (u_int64_t)(999999999));

		_nbPartitions = _options->getInt(STR_SIMKA2_NB_PARTITION);
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
		_nbBankPerDataset = _options->getInt("-nb-dataset");

		//_minKmerShannonIndex = _options->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
		//_minKmerShannonIndex = std::max(_minKmerShannonIndex, 0.0);
		//_minKmerShannonIndex = std::min(_minKmerShannonIndex, 2.0);

		//if(!System::file().doesExist(_inputFilename)){
		//	cerr << "ERROR: Input filename does not exist" << endl;
		//	exit(1);
		//}

		if(!System::file().doesExist(_outputDir)){
			std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
			exit(EXIT_FAILURE);
			/*
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
				std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
				exit(1);
			}*/
		}

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
	}

	void count(){

		IBank* bank = Bank::open(_inputFilename);
		LOCAL(bank);

		SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);
		IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _maxNbReads, _nbBankPerDataset);

		LOCAL(filteredBank);


		Iterator<Sequence>* itSeq = createIterator<Sequence> (
				filteredBank->iterator(),
				filteredBank->estimateNbItems(),
				"Computing minhash sketch and counting"
		);
		LOCAL(itSeq);


		IDispatcher* dispatcher = new SerialDispatcher();
		u_int64_t bloomMemoryBits = _maxMemory * MBYTE * 8;
		Bloom<Type>*  bloomFilter = new BloomCacheCoherent<Type>(bloomMemoryBits, 7);
		//mutex commandMutex;
		//std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> kmerCountSorter;
		//unordered_map<u_int64_t, KmerCountType> kmerCounts;
		SelectKmersCommand<span> command(_kmerSize, _sketchSize, _kmerCountSorter, _kmerCounts, bloomFilter);

		dispatcher->iterate (itSeq, command, 1000);

		delete dispatcher;
		delete bloomFilter;

	}


	void partitionKmerCounts(){

		_partitionWriter = new SimkaPartitionWriter<span>(_outputDir, _nbPartitions);

		size_t sketchSize = _kmerCountSorter.size();

		vector<u_int64_t> _reverseQueue(sketchSize);
		for(size_t i=0; i<sketchSize; i++){
			u_int64_t kmer = _kmerCountSorter.top();
			_reverseQueue[sketchSize-i-1] = kmer;
			_kmerCountSorter.pop();
		}


		for(size_t i=0; i<sketchSize; i++){
			u_int64_t kmer = _reverseQueue[i];
			u_int32_t count = _kmerCounts[kmer];
			//cout << kmer << " " << count << endl;
			//cout << kmerCount._kmer << "  " << kmerCount._count << endl;
			//_partitionWriter->insert(kmerCount._kmer, _datasetIDbin, kmerCount._count);
			_partitionWriter->insert(kmer, _datasetIDbin, count);
			//cout << _minHashKmers[i] << " " << _minHashKmersCounts[i] << endl;
		}


		_partitionWriter->end();

		delete _partitionWriter;
	}


	void saveInfos(){

		u_int64_t nbDistinctKmers = 0;
		u_int64_t nbKmers = 0;
		u_int64_t chord_N2 = 0;
		for(size_t i=0; i<_nbPartitions; i++){
			nbDistinctKmers += _partitionWriter->_nbDistinctKmerPerParts[i];
			nbKmers += _partitionWriter->_nbKmerPerParts[i];
			chord_N2 += _partitionWriter->_chordNiPerParts[i];
		}

        ofstream outputInfoFile((_outputDir + "/merge_info.bin").c_str(), std::ios::binary);

        u_int64_t nbDatasets = 1;
        u_int64_t datasetIdSize = _datasetID.size();
        u_int64_t kmerSize = _kmerSize;

        outputInfoFile.write((char const*)(&nbDatasets), sizeof(nbDatasets));
        SimkaIoUtils::simka2_writeDatasetInfo(outputInfoFile, _datasetID, _nbReads, nbDistinctKmers, nbKmers, chord_N2);
        /*
        outputInfoFile.write((char const*)(&_datasetIDbin), sizeof(_datasetIDbin));
        outputInfoFile.write((char const*)(&datasetIdSize), sizeof(datasetIdSize));
        outputInfoFile.write((char const*)(&_datasetID), _datasetID.size());
        outputInfoFile.write((char const*)(&kmerSize), sizeof(kmerSize));
        outputInfoFile.write((char const*)(&_nbReads), sizeof(_nbReads));
        outputInfoFile.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
        outputInfoFile.write((char const*)(&nbKmers), sizeof(nbKmers));
        outputInfoFile.write((char const*)(&chord_N2), sizeof(chord_N2));
        */

        outputInfoFile.close();
	}


};











class Simka2ComputeKmerSpectrum : public Tool{
public:

	string _execFilename;

	Simka2ComputeKmerSpectrum(string execFilename): Tool ("Simka2-ComputeKmerSpectrum"){
		_execFilename = execFilename;

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

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
	    kmerParser->push_back (new OptionOneParam (STR_SIMKA_SKETCH_SIZE, "number of kmers used to compute distances", true));
	    //kmerParser->push_back(dskParser->getParser (STR_KMER_SIZE));
	    //kmerParser->push_back(new OptionOneParam (STR_KMER_PER_READ.c_str(), "number of selected kmers per read", false, "0"));
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "1"));
	    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "2"));
	    KmerCountType maxAbundance = -1;
	    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "max abundance a kmer can have to be considered", false, Stringify::format("%i", maxAbundance)));

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
	    readParser->push_back (new OptionOneParam ("-nb-dataset", "nb paired datasets", true));

	    //Core parser
	    IOptionsParser* coreParser = new OptionsParser ("core");
	    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
	    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (MB)", false, "8000"));
	    coreParser->push_back (new OptionOneParam (STR_SIMKA2_NB_PARTITION, "nb partitions", true));


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
		if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE))) { p->setDefaultValue ("7"); }
		dskParser->setVisible(false);

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
	    Parameter (IProperties* props, const string& execFilename) : _props(props), _execFilename(execFilename) {}
	    //Simka& _simka;
	    IProperties* _props;
	    string _execFilename;
	};

	template<size_t span> struct Functor {

    	void operator ()  (Parameter p){

    		Simka2ComputeKmerSpectrumAlgorithm<span>* algo = new Simka2ComputeKmerSpectrumAlgorithm<span>(p._props, p._execFilename);
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
    	Simka2ComputeKmerSpectrum(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



