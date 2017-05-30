

#include <gatb/gatb_core.hpp>
#include "ProbabilisticDictionary.hpp"
#include "SimkaMinUtils.hpp"
#include "SimkaMinDistance.hpp"
#include "MurmurHash3.h"
#include <unordered_set>

//SimkaMin n'effectue pas assez de calcul par rapport au temps de decompression des fichiers .gz qui est faites ens érie actuellement dans GATB
//La parallelisation de Simka est donc faites au niveau des jeux de données
//Cependant, la parallelisation au niveau d'un seul jeu de données est implémenter grace a un dispatcher de reads de GATB
//CORE_PER_THREAD indique combien de cores alloué à ce dispatcher
//Le maximum de CPU obtenu par le dispatcher est d'environ 400% sur un jeu .gz, il ne faut donc pas monter cette valeur au dela de 4
//Pas le mettre trop bas tout de même car cette parallelization interne est bonne car elle réduit le nombre de jeu lu en parallele
#define CORE_PER_THREAD 1 //Core for counting k-mer in a given dataset, do not put to high value because decompressing gz is so slow compared to computation

const string STR_SIMKA_SOLIDITY_PER_DATASET = "-solidity-single";
const string STR_SIMKA_MAX_READS = "-max-reads";
const string STR_SIMKA_MIN_READ_SIZE = "-min-read-size";
const string STR_SIMKA_MIN_READ_SHANNON_INDEX = "-read-shannon-index";
const string STR_SIMKA_MIN_KMER_SHANNON_INDEX = "-kmer-shannon-index";
const string STR_KMER_PER_READ = "-kmer-per-read";
const string STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES= "-simple-dist";
const string STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES = "-complex-dist";
const string STR_SIMKA_KEEP_TMP_FILES = "-keep-tmp";
const string STR_SIMKA_COMPUTE_DATA_INFO = "-data-info";

const string STR_SIMKA_SUBSAMPLING_SETUP = "-subsampling-setup";
const string STR_SIMKA_SUBSAMPLING_MAX_READS = "-subsampling-space";
const string STR_SIMKA_SUBSAMPLING_NB_PICKED_READS = "-subsampling-nb-reads";
const string STR_SIMKA_SUBSAMPLING_REFERENCE_DATASET_ID = "-subsampling-ref-id";
const string STR_SIMKA_SUBSAMPLING_KIND = "-subsampling-kind";
const string STR_SIMKA_NB_KMERS_USED = "-nb-kmers";


typedef ProbabilisticDictionary<u_int32_t> ProbabilisticDict;
typedef unordered_map<u_int64_t, u_int32_t> SelectedKmerIndexType;

/*
reprise:
- selectionner les avec l'approche minhash, on se limite donc a un nombre de kmer selectionner faible, max 1M
- selectionner les kmers solide avec un bloom filter (utiliser la meme approche que mash), il va donc probablement falloir paralleliser la selection au niveau des jeux de donnée par rapport à ce filtre
- indexer et comtper les kmers selectionner avec une hasmap (unordered_map), vu qu'on est sur de petit ensemble de kmers
- prendre en compte la mémoire des sketch minHash dans le maxmemory (actuellement tout est donner au bloom filter)
- tester de selectionner les kmer (minhash) que sur une fraction du maxReads. peut etre que tout ce qui importe c'est la fonction de hachage murmurhash...


			reprise: filtrer les kmers avec le _bloomFilter
			virer le quasidict remplacer par une unordered_map
*/



//reprise: verifier que ça marche bien le serial dispatcher par rapport au destructeur du functor...


















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
	mutex& _mutex;

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
	unordered_set<u_int64_t> _kmerCounts;

	std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter>& _kmerCountSorterSynch;
	unordered_set<u_int64_t>& _kmerCountsSynch;

	Bloom<KmerType>* _bloomFilter;
	u_int64_t _nbInsertedKmersInBloom;

	SelectKmersCommand(size_t kmerSize, size_t sketchSize, mutex& mutex, std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter>& kmerCountSorterSynch, unordered_set<u_int64_t>& kmerCountsSynch, Bloom<KmerType>* bloomFilter)
	: _model(kmerSize), _itKmer(_model), _mutex(mutex), _kmerCountSorterSynch(kmerCountSorterSynch), _kmerCountsSynch(kmerCountsSynch), _bloomFilter(bloomFilter)
	{
		_kmerSize = kmerSize;
		_sketchSize = sketchSize;
		_isMaster = true;
		_nbInsertedKmersInBloom = 0;
	}

	SelectKmersCommand(const SelectKmersCommand& copy)
	: _model(copy._kmerSize), _itKmer(_model), _mutex(copy._mutex), _kmerCountSorterSynch(copy._kmerCountSorterSynch), _kmerCountsSynch(copy._kmerCountsSynch), _bloomFilter(copy._bloomFilter)
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



		_mutex.lock();

		cout << "Inserted k-mers: " << _nbInsertedKmersInBloom << endl;
		//cout << "deleteeeeee" << endl;
		//cout << this << endl;

		_kmerCountSorter.pop(); //Discard greater element because queue size is always equal to (_sketchSize + 1) because of an optimization

		size_t sketchSize = _kmerCountSorter.size();
		for(size_t i=0; i<sketchSize; i++){
			u_int64_t kmer = _kmerCountSorter.top();
			_kmerCountSorter.pop();
			mergeSynch(kmer);
			//KmerCount kmerCount = _kmerCountSorter.top();
			//cout << kmerCount._kmer << "  " << kmerCount._count << endl;
			//_partitionWriter->insert(kmerCount._kmer, _datasetIDbin, kmerCount._count);
			//_kmerCountSorter.pop();
			//_partitionWriter->insert(_minHashKmers[i], _datasetIDbin, _minHashKmersCounts[i] );
			//cout << _minHashKmers[i] << " " << _minHashKmersCounts[i] << endl;
		}

		//cout << "deleteeeeee1" << endl;

		_mutex.unlock();
	}

	inline void mergeSynch(u_int64_t kmer){
		if(_kmerCountsSynch.find(kmer) == _kmerCountsSynch.end()){
			if(_kmerCountSorterSynch.size() > _sketchSize){
				if(kmer < _kmerCountSorterSynch.top() ){
					//cout << kmer << "     " << _kmerCounts.size() << endl;
					u_int64_t greaterValue = _kmerCountSorterSynch.top();
					_kmerCountsSynch.erase(greaterValue);
					_kmerCountSorterSynch.pop();
					_kmerCountSorterSynch.push(kmer);
					_kmerCountsSynch.insert(kmer);
				}
				//else{
				//	cout << "\t\tnonon" << endl;
				//}
			}
			else{ //Filling the queue with first elements
				_kmerCountSorterSynch.push(kmer);
				_kmerCountsSynch.insert(kmer);
			}
		}
		//else{
		//	_kmerCountsSynch[kmer] += count;
		//}
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
			MurmurHash3_x64_128 ( (const char*)&kmerValue, sizeof(kmerValue), 100, &kmerHashed);

			//todo: verifier dabord si le kmer peut etre insérer, plus rapide que els accès au table de hachage (bloom et selected)

			if(_kmerCountSorter.size() < _sketchSize){
				if(_bloomFilter->contains(kmer)){
					//Filling the queue with first elements
					if(_kmerCounts.find(kmerHashed) == _kmerCounts.end()){
						_kmerCountSorter.push(kmerHashed);
						_kmerCounts.insert(kmerHashed);
						//cout << _kmerCountSorter.size() << endl;
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
							_kmerCounts.insert(kmerHashed);
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
class CountKmerCommand
{
public:


    typedef typename Kmer<span>::Type  KmerType;
	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;

	size_t _kmerSize;
	ModelCanonical _model;
	ModelCanonicalIterator _itKmer;
	//ProbabilisticDict* _selectedKmersIndex;
	SelectedKmerIndexType& _selectedKmersIndex;

	vector<mutex>* _master_mutex;
	vector<KmerCountType>* _master_kmerCounts;

	vector<mutex>* _mutex;
	vector<KmerCountType>* _kmerCounts;

	bool _isMaster;
	bool _exists;
	size_t _nbMutex;
	//vector<mutex>& _synchMutex;
	//vector<KmerCountType>* _values;
	KmerCountType _valuesMax;


	CountKmerCommand(size_t kmerSize, SelectedKmerIndexType& selectedKmersIndex, size_t nbCores, u_int64_t nbElements)
	: _model(kmerSize), _itKmer(_model), _selectedKmersIndex(selectedKmersIndex)
	{
		_kmerSize = kmerSize;
		_isMaster = true;
		//_master = this;

		cout << nbElements << " " << selectedKmersIndex.size() << endl;
		_nbMutex = nbCores*100;
		_master_mutex = new vector<mutex>(_nbMutex);
		_master_kmerCounts = new vector<KmerCountType>(nbElements, 0);
		_valuesMax = -1; //store the maximum value of ValueType to prevent overflow
	}

	CountKmerCommand(const CountKmerCommand& copy)
	: _model(copy._kmerSize), _itKmer(_model), _selectedKmersIndex(copy._selectedKmersIndex)
	{
		_kmerSize = copy._kmerSize;
		_isMaster = false;

		_valuesMax = -1; //store the maximum value of ValueType to prevent overflow
		_nbMutex = copy._nbMutex;

		_mutex = copy._master_mutex;
		_kmerCounts = copy._master_kmerCounts;

	}

	~CountKmerCommand(){
		if(_isMaster){
			delete _master_mutex;
			delete _master_kmerCounts;
		}
	}


	void operator()(Sequence& sequence){
		//cout << sequence.getIndex() << endl;
		//cout << sequence.toString() << endl;

		_itKmer.setData(sequence.getData());

		for(_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			//cout << _itKmer->value().toString(_kmerSize) << endl;


			KmerType kmer = _itKmer->value();
			u_int64_t kmerValue = kmer.getVal();

			//cout << "1" << endl;
			u_int64_t kmerValueHashed;
			//u_int64_t kmerValue = _itKmer->value().getVal();
			MurmurHash3_x64_128 ( (const char*)&kmerValue, sizeof(kmerValue), 100, &kmerValueHashed);

			//cout << "2" << endl;
			//_mutex->at(0).lock();
			//cout << kmerValueHashed << endl;
			//_mutex->at(0).unlock();

			//cout << "lol" << endl;
			auto element = _selectedKmersIndex.find(kmerValueHashed);

			//_exists = _selectedKmersIndex.find(kmerValueHashed) != _selectedKmersIndex.end();
			//u_int64_t index = _selectedKmersIndex[kmerValueHashed];
			//if(index >= _kmerCounts->size()) cout << "lol" << endl;
			//if(index > 1000){
			//	cout << kmerValue << " " << index  << endl;
			//}
			if(element != _selectedKmersIndex.end()){
				incrementCount(element->second); //todo tester ça
				//_selectedKmersIndex->increment(index);
				//cout << kmerValue << " " << index  << endl;
			}
			//cout << "3" << endl;

		}
	}

	inline void incrementCount(u_int64_t index){
		//cout << index << " " << _kmerCounts->size() << endl;
		_mutex->at(index%_nbMutex).lock();
		u_int64_t newValue = _kmerCounts->at(index)+1;
		//if(index%_nbMutex ==0) cout << _values->at(index) << endl;
		if(newValue < _valuesMax){
			_kmerCounts->at(index) += 1;
		}
		//if(index%_nbMutex ==0)  cout << "\t" << _values->at(index) << endl;
		//this->_values[index].push_back(value);
		_mutex->at(index%_nbMutex).unlock();
	}

};


































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
class SimkaMinAlgorithm : public Algorithm
{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type  KmerType;
	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;
	typedef typename SelectKmersCommand<span>::KmerCountSorter KmerCountSorter;


    u_int64_t _maxMemory;
	size_t _nbCores;
	string _outputDir;
	string _outputDirTemp;
	size_t _nbBanks;
	string _inputFilename;
	size_t _kmerSize;
	int64_t _maxNbReads;
    u_int64_t _nbUsedKmers;
    u_int64_t _nbUsedKmersPerDataset;
    KmerCountType _abundanceMin;
    KmerCountType _abundanceMax;

    size_t _nbCoresPerThread;
	//size_t _minReadSize;
	//double _minReadShannonIndex;
	//double _minKmerShannonIndex;
	//size_t _nbMinimizers;
	//size_t _nbCores;

	//SimkaStatistics* _stats;
	//SimkaDistance* _simkaDistance;

	string _banksInputFilename;
	//vector<string> _tempFilenamesToDelete;
	//IBank* _banks;
	IProperties* _args;

	vector<string> _bankNames;
	vector<size_t> _reorderedBankIds;
	ofstream _outputKmerCountFile;
	//vector<u_int64_t> _nbReadsPerDataset;

	//string _outputFilenameSuffix;

	//u_int64_t _totalKmers;
	vector<size_t> _nbBankPerDataset;
	u_int64_t _bloomFilterMaxMemoryBits;

	//string _largerBankId;
	//bool _computeSimpleDistances;
	//bool _computeComplexDistances;
	//bool _keepTmpFiles;


	SelectedKmerIndexType _selectedKmersIndex;
	//Bloom<Type>*  _bloomFilter;
	SimkaMinDistance<KmerCountType> _distanceManager;

	//Bloom<Type>*  _bloomFilter;
	//Bloom<Type>*  _bloomFilter2;


	std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> _selectedKmerSorter;
	unordered_set<u_int64_t> _selectedKmersSet;


    SimkaMinAlgorithm(IProperties* args):
		Algorithm("simkaMin", -1, args)
	{

    	_args = args;
	}

	void execute(){
		parseArgs();
		createDirs();
		layoutInputFilename();
		checkInputsValidity();

		//computeNbUsedKmers();
		selectKmers();

		//reprise: la stratégie actuelle est trop lente car pas assez e calcul par rapport au decompression gz. abandonner l'idée du dispatcher et paralleliser au niveau des jeux de données
		//string filename = _outputDirTemp + "/selectedKmers.bin";
		//file_binary<u_int64_t> keyFile(filename.c_str());
		//todo load hash table
		//_selectedKmersIndex = new ProbabilisticDict(_nbUsedKmers, keyFile, _nbCores);

		countKmers();
		computeDistance();
	}

	void parseArgs(){
		KmerCountType maxAbundance = -1;
		//_keepTmpFiles = _options->get(STR_SIMKA_KEEP_TMP_FILES);
		_maxMemory = _args->getInt(STR_MAX_MEMORY);
		_nbCores = _args->getInt(STR_NB_CORES);
		_inputFilename = _args->getStr(STR_URI_INPUT);
		_outputDir = _args->get(STR_URI_OUTPUT) ? _args->getStr(STR_URI_OUTPUT) : "./";
		_outputDirTemp = _args->get(STR_URI_OUTPUT_TMP) ? _args->getStr(STR_URI_OUTPUT_TMP) : "./";
		_kmerSize = _args->getInt(STR_KMER_SIZE);
		_nbUsedKmers = _args->getInt(STR_SIMKA_NB_KMERS_USED);
		_abundanceMin = _args->getInt(STR_KMER_ABUNDANCE_MIN);
		_abundanceMax = min((u_int64_t)_args->getInt(STR_KMER_ABUNDANCE_MAX), (u_int64_t)maxAbundance);

		//cout << _options->getInt(STR_KMER_ABUNDANCE_MAX) << endl;
		//cout << _abundanceThreshold.second << endl;
		//_soliditySingle = _options->get(STR_SIMKA_SOLIDITY_PER_DATASET);
		//_nbMinimizers = _options->getInt(STR_KMER_PER_READ);
		//_maxDisk = getInput()->getInt(STR_MAX_DISK);

		//read filter
		_maxNbReads = _args->getInt(STR_SIMKA_MAX_READS);
		//_minReadSize = _options->getInt(STR_SIMKA_MIN_READ_SIZE);
		//_minReadShannonIndex = _options->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
		//_minReadShannonIndex = std::max(_minReadShannonIndex, 0.0);
		//_minReadShannonIndex = std::min(_minReadShannonIndex, 2.0);

		//_minKmerShannonIndex = _options->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
		//_minKmerShannonIndex = std::max(_minKmerShannonIndex, 0.0);
		//_minKmerShannonIndex = std::min(_minKmerShannonIndex, 2.0);

		if(!System::file().doesExist(_inputFilename)){
			cerr << "ERROR: Input filename does not exist" << endl;
			exit(1);
		}


		_nbCoresPerThread = min((u_int64_t)_nbCores, (u_int64_t)CORE_PER_THREAD);

	}

	void createDirs(){

		if(!System::file().doesExist(_outputDir)){
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
		        std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
		        exit(1);
			}
		}

		if(!System::file().doesExist(_outputDirTemp)){
			int ok = System::file().mkdir(_outputDirTemp, -1);
			if(ok != 0){
		        std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
		        exit(1);
			}
		}

		_outputDirTemp = System::file().getRealPath(_outputDirTemp);
		_outputDirTemp += "/simka_output_temp/";
		System::file().mkdir(_outputDirTemp, -1);

		_args->setStr(STR_URI_OUTPUT_TMP, _outputDirTemp);
		System::file().mkdir(_outputDirTemp + "/input/", -1);

	}


	void layoutInputFilename(){

		string inputDir = _outputDirTemp + "/input/";
		ifstream inputFile(_inputFilename.c_str());

		_banksInputFilename =  inputDir + "__input_simka__"; //_inputFilename + "_dsk_dataset_temp__";
		IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		string line;
		string linePart;
		vector<string> lineIdDatasets;
		vector<string> linepartPairedDatasets;
		vector<string> linepartDatasets;

		string bankFileContents = "";

		u_int64_t lineIndex = 0;

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
			_nbBankPerDataset.push_back(linepartPairedDatasets.size());

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

			bankFileContents += inputDir + "/" + bankId + "\n";
			lineIndex += 1;

			_bankNames.push_back(bankId);


		}


		inputFile.close();

		bankFileContents.erase(bankFileContents.size()-1);
		bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);
		bankFile->flush();
		delete bankFile;

		_nbBanks = _bankNames.size();

		if(_args->getInt(STR_VERBOSE) != 0){
			cout << "Nb input datasets: " << _nbBanks << endl;
			cout << endl;
		}

	}


	void checkInputsValidity(){

		bool error = false;
		string inputDir = _outputDirTemp + "/input/";

		for (size_t i=0; i<_nbBanks; i++){

			try{
				IBank* bank = Bank::open(inputDir + _bankNames[i]);
				LOCAL(bank);
			}
			catch (Exception& e){
				cerr << "ERROR: Can't open dataset: " << _bankNames[i] << endl;
				error = true;
			}

		}

		if(error) exit(1);
	}

	/*
	void computeNbUsedKmers(){
		cout << endl << endl;
		u_int64_t bloomFilterMaxMemoryByte = _maxMemory*GBYTE* (2.0/3.0);
		u_int64_t bloomFilterMaxMemoryBits = bloomFilterMaxMemoryByte*8;
		//float bytePerKmers = (float)nbBitPerKmers/8.0;
		u_int64_t maxUsableKmers = (((long double)bloomFilterMaxMemoryBits)/_nbBitPerKmers);

		cout << "Nb kmers wanted  : " << _nbUsedKmers << endl;
		cout << "Max usable kmers : " << maxUsableKmers << endl;

		if(_nbUsedKmers == 0) _nbUsedKmers = maxUsableKmers;

		_nbUsedKmers = min(maxUsableKmers, _nbUsedKmers);
		cout << "Nb used kmers    : " << _nbUsedKmers << endl;
		_nbUsedKmersPerDataset = _nbUsedKmers / _nbBanks;
		cout << "Nb used kmers per dataset : " << _nbUsedKmersPerDataset << endl;
		cout << endl << endl;

		_bloomFilterMaxMemoryBits = _nbUsedKmers * _nbBitPerKmers;
		cout << "Bloom filter size (bits) : " << _bloomFilterMaxMemoryBits << endl;
	}


	void selectKmers(){

		string filename = _outputDirTemp + "/selectedKmers.bin";
		ofstream selectKmersFile(filename.c_str(), ios::binary);

		//todo trouver la bonne taille de filtre de bloom en fonction des params
		u_int64_t mainBloomFilterMemoryBits = _maxMemory*GBYTE*0.5*8;
		_bloomFilter = new BloomCacheCoherent<Type> (mainBloomFilterMemoryBits, 7);
		u_int64_t subBloomFilterMemoryBits = _maxMemory*GBYTE*0.5*8;

		IBank* bank = Bank::open(_banksInputFilename);
		LOCAL(bank);


		ModelCanonical model(_kmerSize);
		ModelCanonicalIterator itKmer(model);

		//SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);
		//IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _maxNbReads, _nbBankPerDataset);

		//LOCAL(filteredBank);

		Iterator<Sequence>* it = bank->iterator();
		//Iterator<Sequence>* itSeq = createIterator<Sequence> (
		//		bank->iterator(),
		//		bank->estimateNbItems(),
		//		"Selecting and indexing kmers"
		//);
		//LOCAL(itSeq);

	    std::vector<Iterator<Sequence>*> itBanks =  it->getComposition();
	    u_int64_t nbKmerInsertedTotal = 0;
	    size_t skip = 0;

	    for (size_t i=0; i<itBanks.size(); i++)
	    {

			_bloomFilter2 = new BloomCacheCoherent<Type> (subBloomFilterMemoryBits, 7);

	    	Iterator<Sequence>* itSeq = itBanks[i];
	    	u_int64_t nbKmerInserted = 0;
	    	bool isDone = false;

			for(itSeq->first(); !itSeq->isDone(); itSeq->next()){
				Sequence& sequence = itSeq->item();

				itKmer.setData(sequence.getData());

				for(itKmer.first(); !itKmer.isDone(); itKmer.next()){

					if(skip > 0){
						skip -= 1;
						continue;
					}

					Type kmer = itKmer->value();
					if(kmer.getVal() == 0) continue; //todo can be optimized maybe ?
					if(getShannonIndex(kmer) < 1.83) continue;

					if(_bloomFilter->contains(kmer)){
						skip = _kmerSize;
					}
					else if(_bloomFilter2->contains(kmer)){
						u_int64_t kmerValue = kmer.getVal();
						//cout << nbKmerInsertedTotal << ": " << kmer.toString(_kmerSize) << ": " << getShannonIndex(kmer) << endl;
						selectKmersFile.write((const char*)&kmerValue, sizeof(kmerValue));
						_bloomFilter->insert(kmer);
						nbKmerInserted += 1;
						nbKmerInsertedTotal += 1;
						if(nbKmerInserted >= _nbUsedKmersPerDataset){
							isDone = true;
							break;
						}
						skip = _kmerSize;
					}
					else{
						_bloomFilter2->insert(kmer);
					}
				}

				if(isDone){
					cout << nbKmerInserted << endl;
					break;
				}

			}

			delete _bloomFilter2; //todo could be better if we could clear the bit buffer
			itBanks[i]->finalize();
	    }

	    _nbUsedKmers = nbKmerInsertedTotal;
	    cout << "Final number of selected kmers: " << _nbUsedKmers << endl;
	    selectKmersFile.close();
	    delete _bloomFilter;
	}*/

	vector<thread*> _threads;
	size_t _maxRunningThreads;
	vector<size_t> _runningThreadIds;
	size_t _nbRunningThreads;
	vector<size_t> _finishedThreads;
	mutex countKmersMutex;

	/*
	void obtainThreadId(){
		for(size_t i=0; i<_isThreadRunning.size(); i++){
			if(!_isThreadRunning[i]){
				_isThreadRunning[i] = true;
				return i;
			}
		}
	}*/

	u_int64_t _availableMemoryBitsPerThread;

	void selectKmers(){

		cout << endl << endl;
		cout << "Selecting kmers..." << endl;

		_maxRunningThreads = _nbCores;

		u_int64_t nbBitPerKmers = 12;
		u_int64_t availableMemoryBits = _maxMemory*MBYTE*(u_int64_t)8;
		_availableMemoryBitsPerThread = availableMemoryBits / _maxRunningThreads;
		cout << "Memory allowed to bloom filter (bits): " << _availableMemoryBitsPerThread << endl;
		cout << "bits per k-mer: " << nbBitPerKmers << endl;
		cout << "Max kmer inserted in bloom filter: " << (_availableMemoryBitsPerThread / nbBitPerKmers) << endl;


		//cout << _maxRunningThreads << endl;

		size_t threadId = 0;
		//vector<thread> threads; //(_nbCores);
		//_isThreadRunning = vector<bool>(_nbCores);
		_nbRunningThreads = 0;

		for (size_t i=0; i<_nbBanks; i++){
			thread* t = new thread(&SimkaMinAlgorithm::selectKmersOfDataset, this,  i, threadId);
			_threads.push_back(t);
			_runningThreadIds.push_back(threadId);
			threadId += 1;
			_nbRunningThreads += 1;
			//_isThreadRunning[threadId] = true;
			//_nbRunningThreads[i] += 1;

			if(_nbRunningThreads >= _maxRunningThreads){
				waitThreads();
			}


		}

		joinThreads();

		//string filename = _outputDirTemp + "/selectedKmers.bin";
		//ofstream selectKmersFile(filename.c_str(), ios::binary);

		cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;
		_selectedKmerSorter.pop(); //there is always one extra element because of a >= optimization...
		cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;
		u_int64_t size = _selectedKmerSorter.size();
		for(size_t i=0; i<size; i++){
			u_int64_t kmerValue = _selectedKmerSorter.top();
			_selectedKmerSorter.pop();
			_selectedKmersIndex[kmerValue] = i;
			//selectKmersFile.write((const char*)&kmerValue, sizeof(kmerValue)); //todo mphf loading can be done in memory, not required to write all selected kmers on disk
		}
		_nbUsedKmers = _selectedKmersIndex.size() ;
		cout << _selectedKmerSorter.size() << " " << _nbUsedKmers << endl;

		//selectKmersFile.close();
	}


	void selectKmersOfDataset(size_t datasetId, size_t threadId){

		countKmersMutex.lock();
		IBank* bank = Bank::open(_outputDirTemp + "/input/" + _bankNames[datasetId]);
		LOCAL(bank);
		Iterator<Sequence>* itSeq = bank->iterator();
		LOCAL(itSeq);
		Iterator<Sequence>* itSeqSimka = new SimkaInputIterator<Sequence> (itSeq, _nbBankPerDataset[datasetId], _maxNbReads);
		LOCAL(itSeqSimka);

		Bloom<Type>*  bloomFilter = new BloomCacheCoherent<Type> (_availableMemoryBitsPerThread, 7);

		IDispatcher* dispatcher = new SerialDispatcher();
		mutex commandMutex;
		std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> kmerCountSorter;
		unordered_set<u_int64_t> kmerCounts;
		SelectKmersCommand<span> command(_kmerSize, _nbUsedKmers, commandMutex, kmerCountSorter, kmerCounts, bloomFilter);
		//CountKmerCommand<span> command(_kmerSize, _selectedKmersIndex, _nbCoresPerThread, _nbUsedKmers);



		countKmersMutex.unlock();




		dispatcher->iterate(itSeqSimka, command, 1000);



		countKmersMutex.lock();

		delete bloomFilter;

		//kmerCountSorter.pop(); //Discard greater element because queue size is always equal to (_sketchSize + 1) because of an optimization

		u_int64_t size = kmerCountSorter.size();
		for(size_t i=0; i<size; i++){
			u_int64_t kmer = kmerCountSorter.top();
			kmerCountSorter.pop();

			if(_selectedKmersSet.find(kmer) == _selectedKmersSet.end()){
				if(_selectedKmerSorter.size() > _nbUsedKmers){
					if(kmer < _selectedKmerSorter.top() ){
						//cout << kmer << "     " << _kmerCounts.size() << endl;
						u_int64_t greaterValue = _selectedKmerSorter.top();
						_selectedKmersSet.erase(greaterValue);
						_selectedKmerSorter.pop();
						_selectedKmerSorter.push(kmer);
						_selectedKmersSet.insert(kmer);
					}
					//else{
					//	cout << "\t\tnonon" << endl;
					//}
				}
				else{ //Filling the queue with first elements
					_selectedKmerSorter.push(kmer);
					_selectedKmersSet.insert(kmer);
				}
			}
			//else{
			//	_kmerCountsSynch[kmer] += count;
			//}
			//mergeSynch(kmer, _kmerCounts[kmer]);
			//KmerCount kmerCount = _kmerCountSorter.top();
			//cout << kmerCount._kmer << "  " << kmerCount._count << endl;
			//_partitionWriter->insert(kmerCount._kmer, _datasetIDbin, kmerCount._count);
			//_kmerCountSorter.pop();
			//_partitionWriter->insert(_minHashKmers[i], _datasetIDbin, _minHashKmersCounts[i] );
			//cout << _minHashKmers[i] << " " << _minHashKmersCounts[i] << endl;
		}



		//cout << "done" << endl;
		/*
		//_reorderedBankIds.push_back(datasetId);

		u_int64_t nbKmersWithCounts = 0;
		vector<KmerCountType>* counts = command._master_kmerCounts;
		for(size_t i=0; i<counts->size(); i++){

			KmerCountType kmerCount = counts->at(i);

			if(kmerCount > 0){
				nbKmersWithCounts += 1;
			}
			_outputKmerCountFile.seekp(i*_nbBanks*sizeof(KmerCountType) + sizeof(KmerCountType)*datasetId);
			_outputKmerCountFile.write((const char*)&kmerCount, sizeof(kmerCount));
		}
		cout << "Fill rate of kmer counts: " << ((double)nbKmersWithCounts / (double)counts->size()) << endl;
		//cout << "\t thread " <<  threadId << " End" << endl;*/
		_finishedThreads.push_back(threadId);
		countKmersMutex.unlock();
	}





	void countKmers(){

		cout << endl << endl;
		cout << "Counting selected kmers..." << endl;


		string filename = _outputDirTemp + "/kmerCounts.bin";
		_outputKmerCountFile.open(filename.c_str(), ios::binary);

		u_int64_t nbCounts = _nbUsedKmers*_nbBanks;
		u_int16_t count0 = 0;
		for(size_t i=0; i<nbCounts; i++){
			_outputKmerCountFile.write((const char*)&count0, sizeof(count0));
		}
		_outputKmerCountFile.flush();

		//IBank* bank = Bank::open(_banksInputFilename);
		//LOCAL(bank);
		//Iterator<Sequence>* it = bank->iterator();
		//std::vector<Iterator<Sequence>*> itBanks =  it->getComposition();

		//cout << _nbCores << endl;
		_maxRunningThreads = _nbCores; // / _nbCoresPerThread;
		//_maxRunningThreads = 1;
		//cout << _maxRunningThreads << endl;

		size_t threadId = 0;
		//vector<thread> threads; //(_nbCores);
		//_isThreadRunning = vector<bool>(_nbCores);
		_nbRunningThreads = 0;

		for (size_t i=0; i<_nbBanks; i++){
			//cout << "processing bank: " << i << endl;

			//LOCAL(itSeqSimka);

			//size_t threadId = obtainThreadId();

			thread* t = new thread(&SimkaMinAlgorithm::countKmersOfDataset, this,  i, threadId);
			_threads.push_back(t);
			_runningThreadIds.push_back(threadId);
			threadId += 1;
			_nbRunningThreads += 1;
			//_isThreadRunning[threadId] = true;
			//_nbRunningThreads[i] += 1;

			if(_nbRunningThreads >= _maxRunningThreads){
				waitThreads();
			}



			//---------------------------
			//todo !!!!!!!!!!!!!!!!!!! delete itSeqSimka, thread, itBanks[i]->finalize();
			//---------------------------
		}

		joinThreads();

		_outputKmerCountFile.close();
		//delete _selectedKmersIndex;
		/*
		cout << endl << endl;
		cout << "Start counting selected kmers" << endl;




		for (size_t i=0; i<itBanks.size(); i++)
		{

			Iterator<Sequence>* itSeq = itBanks[i];
			Iterator<Sequence>* itSeqSimka = new SimkaInputIterator<Sequence> (itSeq, _nbBankPerDataset[i], _maxNbReads);
			LOCAL(itSeqSimka);
			//IBank* filteredBank = new SimkaInputIterator<SimkaSequenceFilter>(bank, sequenceFilter, _maxNbReads, _nbBankPerDataset);

			getDispatcher()->iterate(itSeqSimka, CountKmerCommand<span>(_kmerSize, _selectedKmersIndex), 1000);


			itBanks[i]->finalize();
		}*/

	}




	void countKmersOfDataset(size_t datasetId, size_t threadId){

		countKmersMutex.lock();
		//cout << "starting thread: " << threadId << endl;
		IBank* bank = Bank::open(_outputDirTemp + "/input/" + _bankNames[datasetId]);
		LOCAL(bank);
		Iterator<Sequence>* itSeq = bank->iterator();
		LOCAL(itSeq);
		Iterator<Sequence>* itSeqSimka = new SimkaInputIterator<Sequence> (itSeq, _nbBankPerDataset[datasetId], _maxNbReads);
		LOCAL(itSeqSimka);

		//_nbCoresPerThread = 1;
		IDispatcher* dispatcher = new SerialDispatcher();
		//IDispatcher* dispatcher = new Dispatcher (_nbCoresPerThread);
		CountKmerCommand<span> command(_kmerSize, _selectedKmersIndex, _nbCoresPerThread, _nbUsedKmers);

		cout << "starting thread: " << threadId << endl;
		countKmersMutex.unlock();

		dispatcher->iterate(itSeqSimka, command, 1000);

		countKmersMutex.lock();
		cout << "done: " << threadId << endl;
		//_reorderedBankIds.push_back(datasetId);

		u_int64_t nbKmersWithCounts = 0;
		vector<KmerCountType>* counts = command._master_kmerCounts;
		for(size_t i=0; i<counts->size(); i++){

			KmerCountType kmerCount = counts->at(i);

			if(kmerCount > 0){
				nbKmersWithCounts += 1;
				_outputKmerCountFile.seekp(i*_nbBanks*sizeof(KmerCountType) + sizeof(KmerCountType)*datasetId);
				_outputKmerCountFile.write((const char*)&kmerCount, sizeof(kmerCount));
			}
		}
		cout << "Fill rate of kmer counts: " << ((double)nbKmersWithCounts / (double)counts->size()) << endl;
		//cout << "\t thread " <<  threadId << " End" << endl;
		_finishedThreads.push_back(threadId);
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

	void computeDistance(){

		cout << endl << endl;
		cout << "Computing distances" << endl;

		_distanceManager = SimkaMinDistance<KmerCountType>(_nbBanks, _abundanceMin, _abundanceMax);
		u_int64_t nbDistinctKmers = 0;

		vector<KmerCountType> counts(_nbBanks, 0);
		KmerCountType count;
		string filename = _outputDirTemp + "/kmerCounts.bin";
		ifstream kmerCountFile(filename.c_str(), ios::binary);
		size_t datasetId = 0;

		while(!kmerCountFile.eof()){
			kmerCountFile.read((char*)&count, sizeof(count));
			counts[datasetId] = count;

			datasetId += 1;
			if(datasetId >= _nbBanks){

				cout << nbDistinctKmers << ": ";
				for(size_t i=0; i<counts.size(); i++){
					cout << counts[i] << " ";
				}
				cout << endl;

				_distanceManager.processAbundanceVector(counts);
				std::fill(counts.begin(), counts.end(), 0);
				datasetId = 0;
				nbDistinctKmers += 1;
			}

		}

		//cout << endl << endl;
		//cout << "Outputting distances" << endl;
		_distanceManager.computeDistanceMatrix(_outputDirTemp, "mat_abundance_braycurtis");
		_distanceManager.writeMatrixASCII(_outputDir, _outputDirTemp, "mat_abundance_braycurtis", _bankNames);

		cout << endl << endl;
		cout << "Result dir: " << _outputDir << endl;
	}

	double getShannonIndex(const Type&  kmer){
		float index = 0;
		//float freq [5];

		vector<float> _freqs(4, 0);

		//char* seqStr = seq.getDataBuffer();

	    for (size_t i=0; i<_kmerSize; i++){
	    	_freqs[kmer[i]] += 1.0;
	    	//seq[sizeKmer-i-1] = bin2NT [(*this)[i]];
	    }

		// Frequency of each letter (A, C, G, T or N)
		//for(size_t i=0; i < seq.size(); i++)
		//	_freqs[nt2binTab[(unsigned char)seq[i]]] += 1.0;

		// Shannon index calculation
		for (size_t i=0; i<_freqs.size(); i++){
			_freqs[i] /= (float) _kmerSize;
			if (_freqs[i] != 0)
				index += _freqs[i] * log (_freqs[i]) / log(2);
		}
		return abs(index);
	}
};










































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
class SimkaMin : public Tool{

public:


	SimkaMin(): Tool ("SimkaMin"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");


	    //IOptionsParser* parser = parent; //new OptionsParser ("Simka");

		//Main parser
	    //parser->push_front (new OptionNoParam (STR_SIMKA_COMPUTE_DATA_INFO, "compute (and display) information before running Simka, such as the number of reads per dataset", false));
	    parser->push_front (new OptionNoParam (STR_SIMKA_KEEP_TMP_FILES, "keep temporary files", false));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for result files (distance matrices)", false, "./simka_results"));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input file of datasets (see README or examples)", true));


	    //parser->push_back (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
		//IOptionsParser* parser = getParser();
		//IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();
		//parser->push_back(dskParser);
		//dskParser->setVisible(false);
		//cout << parser->getParser(STR_NB_CORES) << endl;
		parser->getParser(STR_NB_CORES)->setVisible(false);

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
	    kmerParser->push_back (new OptionOneParam (STR_SIMKA_NB_KMERS_USED, "nb kmer used to compute distances", false, "0"));
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
	    readParser->push_back (new OptionOneParam (STR_SIMKA_MAX_READS.c_str(), "maximum number of reads per sample to process. Can be -1: use all reads. Can be 0: estimate it", false, "-1" ));
	    //readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE.c_str(), "minimal size a read should have to be kept", false, "0" ));
	    //readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(), "minimal Shannon index a read should have to be kept. Float in [0,2]", false, "0" ));

	    //Core parser
	    IOptionsParser* coreParser = new OptionsParser ("core");
	    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
	    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (MB)", false, "8000"));
	    //coreParser->push_back(dskParser->getParser ());
	    //coreParser->push_back(dskParser->getParser (STR_MAX_DISK));

	    //Subsampling parser
	    //IOptionsParser* subsamplingParser = new OptionsParser ("subsampling");
	    //subsamplingParser->push_back(new OptionNoParam(STR_SIMKA_SUBSAMPLING_SETUP, "provide information for correctly subsampling reads", false));
	    //subsamplingParser->push_back(new OptionOneParam(STR_SIMKA_SUBSAMPLING_MAX_READS, "max size of the subsample space (in reads) (provided by Simka with option: " + string(STR_SIMKA_SUBSAMPLING_SETUP) + ")", false));
	    //subsamplingParser->push_back(new OptionOneParam(STR_SIMKA_SUBSAMPLING_REFERENCE_DATASET_ID, "reference dataset ID for subsampling (provided by Simka with option: " + string(STR_SIMKA_SUBSAMPLING_SETUP) + ")", false));
	    //subsamplingParser->push_back(new OptionOneParam(STR_SIMKA_SUBSAMPLING_NB_PICKED_READS, "number of reads to pick randomly (must be smaller than " + string(STR_SIMKA_SUBSAMPLING_MAX_READS) + ")", false));
	    //subsamplingParser->push_back(new OptionOneParam(STR_SIMKA_SUBSAMPLING_KIND, "0: resampling with replacement, 1: resampling without replacement", false, "0"));
	    //subsamplingParser->push_back(new OptionOneParam(STR_SIMKA_SUBSAMPLING_NB_PICKED_KMERS, "number of kmers to pick randomly", false, "0"));
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
		//parser->push_back(subsamplingParser);
		parser->push_back(coreParser);
		//parser->push_back(distanceParser);


		IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();

	    if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE)))  {  p->setDefaultValue ("7"); }
		parser->push_back(dskParser);
		dskParser->setVisible(false);
	    if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

	}


	struct Parameter
	{
	    Parameter (IProperties* props) : _props(props){}

	    IProperties* _props;
	};

	template<size_t span> struct Functor {

    	void operator ()  (Parameter p){

    		SimkaMinAlgorithm<span>* algo = new SimkaMinAlgorithm<span>(p._props);
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









int main (int argc, char* argv[])
{
    try
    {
    	SimkaMin().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



