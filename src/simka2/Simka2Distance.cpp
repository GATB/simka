/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
//#include "../core/SimkaUtils.hpp"
#include "Simka2Utils.hpp"
#include "../core/SimkaUtils.hpp"
#include "Simka2Database.hpp"
//#include "../minikc/MiniKC.hpp"
//#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
//#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"






template<size_t span>
class DatasetMergerDistance : public DiskBasedMergeSort<span>
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    //typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;
	/*
    struct kxp{
		Type _type;
    	StorageItKmerCount<span>* _it;

    	KmerCount_It(){

    	}

    	KmerCount_It(Count& count, StorageItKmerCount<span>* it){
    		_count = count;
    		_it = it;
    	}
    };



	struct kxpcomp { bool operator() (kxp l,kxp r) { return (get<0>(r) < get<0>(l)); } } ;
	 */

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

	size_t _nbBanks;
	//string _outputDir;
	//string _outputFilename;
	//vector<string>& _datasetToMergeDirs;
	//size_t _partitionId;
	Bag<Kmer_BankId_Count>* _outputGzFile;
	Bag<Kmer_BankId_Count>* _cachedBag;
    CountVector _abundancePerBank;
    SimkaCounterBuilderMerge* _solidCounter;
    Type _lastKmer;
    bool _isInit;
	size_t _nbBankThatHaveKmer;
	SimkaStatistics* _stats;
	SimkaCountProcessorSimple<span>* _processor;
	//bool _computeComplexDistances;
	//size_t _nbBanks;
	//vector<string> _currentDatasetIds;
	u_int64_t _nbDistinctKmers;
	u_int64_t _nbDistinctKmersPerPartition;
	u_int64_t _processedKmerRangeMin;
	u_int64_t _processedKmerRangeMax;
	u_int64_t _partitionRangeMin;
	u_int64_t _partitionRangeMax;
	u_int64_t _lol;
	size_t _currentPartitionId;
	//size_t _nbPartitionToProcess;
	Simka2Database& _database;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	size_t _nbNewBanks;
	size_t _kmerSize;
	size_t _sketchSize;
	vector<string>& _kmerSpectrumDirs;
	string& _dirMatrixParts;

	DatasetMergerDistance(size_t nbBanks, size_t partitionId, size_t nbPartitionToProcess, vector<string>& datasetToMergeDirs, map<string, u_int64_t>& idToOrder, u_int64_t nbDistinctKmersPerPartition,
			Simka2Database& database, bool computeSimpleDistances, bool computeComplexDistances, size_t nbNewBanks, size_t kmerSize, size_t sketchSize, vector<string>& kmerSpectrumDirs, string& dirMatrixParts):
		DiskBasedMergeSort<span>(0, datasetToMergeDirs, idToOrder, true), _database(database), _computeSimpleDistances(computeSimpleDistances), _computeComplexDistances(computeComplexDistances), _nbNewBanks(nbNewBanks), _kmerSize(kmerSize), _sketchSize(sketchSize), _kmerSpectrumDirs(kmerSpectrumDirs), _dirMatrixParts(dirMatrixParts)
    {


		_nbBanks = nbBanks;
		_abundancePerBank.resize(_nbBanks, 0);
		_solidCounter = new SimkaCounterBuilderMerge(_abundancePerBank);
		_isInit = false;
		_nbBankThatHaveKmer = 0;
		//_stats = stats;
		//_processor = processor;
		_nbDistinctKmers = 0;
		_lol = 0;
		_nbDistinctKmersPerPartition = nbDistinctKmersPerPartition;
		_currentPartitionId = partitionId*nbPartitionToProcess;
		//_nbPartitionToProcess = nbPartitionToProcess;

		cout << partitionId << " " <<  nbPartitionToProcess << " " <<  _nbDistinctKmersPerPartition << endl;
		_partitionRangeMin = partitionId * nbPartitionToProcess;
		_partitionRangeMax = (partitionId+1) * nbPartitionToProcess;
		_processedKmerRangeMin = _partitionRangeMin * _nbDistinctKmersPerPartition;
		_processedKmerRangeMax = _partitionRangeMax * _nbDistinctKmersPerPartition;
		//_processedKmerRangeMax = _processedKmerRangeMin + _nbDistinctKmersPerPartition;
		cout << partitionId << ": " << _processedKmerRangeMin << " " << _processedKmerRangeMax << endl;

		initPartition();
    }


	void process(Type& kmer, u_int64_t bankId, u_int64_t abundance){

		//return;
    	//cout << kmer.toString(31) << " " << bankId <<  " "<< abundance << endl;
		if(_isInit){
			if(kmer == _lastKmer){
				_solidCounter->increase(bankId, abundance);
				_nbBankThatHaveKmer += 1;
			}
			else{
				insert(_lastKmer, _abundancePerBank, _nbBankThatHaveKmer);
				_solidCounter->init(bankId, abundance);
				_lastKmer = kmer;
				_nbBankThatHaveKmer = 1;
			}
		}
		else{
			_isInit = true;
			_lastKmer = kmer;
			_solidCounter->init(bankId, abundance);
	        _nbBankThatHaveKmer = 1;
		}
		//cout << kmer.toString(31) << " " << bankId << " " << abundance << endl;
		//_cachedBag->insert(Kmer_BankId_Count(kmer, bankId, abundance));
	}

	void insert(const Type& kmer, const CountVector& counts, size_t nbBankThatHaveKmer){

		//cout << _processedKmerRangeMin << " " << _processedKmerRangeMax << endl;
		/*
		cout << kmer.toString(31) << "     ";
		for(size_t i=0; i<counts.size(); i++){
			cout << counts[i] << " ";
		}
		cout << endl;
		cout << "-----" << endl;
*/

		//_stats->_nbDistinctKmers += 1;

		//if(_nbBankThatHaveKmer > 1){
		//	_stats->_nbSharedKmers += 1;
		//}

		//if()
		//if(_nbDistinctKmers >= _processedKmerRangeMin && _nbDistinctKmers < _processedKmerRangeMax){
		//	_processor->process(this->_partitionId, kmer, counts);
		//}
		if(_nbDistinctKmers >= _processedKmerRangeMin){
			if(_nbDistinctKmers < _processedKmerRangeMax){
				_processor->process(this->_partitionId, kmer, counts);
				_lol += 1;
			}
			else{
				cout << "------------ DOOOOONE" << endl;
				this->_isDone = true;
				return;
			}
			if(_lol >= _nbDistinctKmersPerPartition){
			//if(_nbDistinctKmers % _nbDistinctKmersPerPartition == 0){
				endPartition();
				initPartition();
			}
		}

		_nbDistinctKmers += 1;
		/*
		if(_stats->_computeComplexDistances || _nbBankThatHaveKmer > 1){

			if(_nbBankThatHaveKmer > 1){
				_stats->_nbSharedKmers += 1;
			}

			_processor->process(this->_partitionId, kmer, counts);

		}*/
	}


	void end(){
		cout << "\tthe last of us: " << _currentPartitionId << endl;
		//cout << "ENNNNNNNNNNNNNNNNNNNNNNNND" << endl;
		insert(_lastKmer, _abundancePerBank, _nbBankThatHaveKmer);
		delete _solidCounter;

		endPartition();
		//cout << "ENNNNNNNNNNNNNNNNNNNNNNNND2" << endl;
    }

	void initPartition(){
		_lol = 0;
		//if(_currentPartitionId >= _partitionRangeMin && _currentPartitionId < _partitionRangeMax){
			cout << "Init partition: " << _currentPartitionId << endl;
			_stats = new SimkaStatistics(_nbBanks, _nbNewBanks, _computeSimpleDistances, _computeComplexDistances, _kmerSize, _sketchSize);
			simka2_loadStatInfos(_database._dir, _database._uniqKmerSpectrumDirs, _database._entries, _kmerSpectrumDirs, _stats, _database._entriesInfos);
			_processor = new SimkaCountProcessorSimple<span> (_stats, _nbBanks, _nbNewBanks, _kmerSize, 0);
			//}
	}

	void endPartition(){
		//cout << "nb processed kmers: " << _lol << endl;
		//if(_currentPartitionId >= _partitionRangeMin && _currentPartitionId < _partitionRangeMax){
			cout << "End partition: " << _currentPartitionId << endl;
			saveStats();

			delete _stats;
			delete _processor;
			_currentPartitionId += 1;
			//}
	}

	void saveStats(){
		string filename = _dirMatrixParts + "/" + Stringify::format("%i", _currentPartitionId) + ".gz";
		_stats->save(filename);
	}

};



















/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span>
class DistinctMergedKmerCounter : public DiskBasedMergeSort<span>
{

public:

	typedef typename Kmer<span>::Type                                       Type;


	u_int64_t _nbDistinctMergedKmers;

    Type _lastKmer;
    bool _isInit;


	DistinctMergedKmerCounter(size_t partitionId, vector<string>& datasetToMergeDirs, map<string, u_int64_t>& idToOrder):
		DiskBasedMergeSort<span>(partitionId, datasetToMergeDirs, idToOrder, false)
    {
		//cout << partitionId << endl;
		//_lastKmer = Type(0);
		_nbDistinctMergedKmers = 0;
		_isInit = false;
    }

	void process(Type& kmer, u_int64_t bankId, u_int64_t abundance){

		if(_isInit){
			if(kmer == _lastKmer){
			}
			else{
				_nbDistinctMergedKmers += 1;
				_lastKmer = kmer;
			}
		}
		else{
			_isInit = true;
			_lastKmer = kmer;
		}
	}



	void end(){
		_nbDistinctMergedKmers += 1;
    }




};







/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2DistanceAlgorithm : public Algorithm
{

public:

	string _databaseDir;
	//string _inputFilenameUnknown;
	string _outputDir;
	size_t _nbCores;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	size_t _kmerSize;
	size_t _partitionId;
	size_t _nbParitionToProcess;
	size_t _maxDatasets;
	size_t _sketchSize;
	size_t _distanceType;
	u_int64_t _nbDistinctMergedKmers;

	Simka2Database _database;
	vector<string> _allIds;
	vector<string> _allKmerSpectrumDirs;
	vector<string> _newIds;

	vector<string> _kmerSpectrumDirs;

	SimkaStatistics* _stats;
	SimkaCountProcessorSimple<span>* _processor;

	size_t _nbBanks;
	size_t _nbNewBanks;

	string _dirMatrixParts;
	map<string, u_int64_t> _idToOrder;
	//string _dirMatrixMatrixBinary;

	Simka2DistanceAlgorithm(IProperties* options):
		Algorithm("simka", -1, options)
	{

	}

	void execute(){
		parseArgs();
		createDatabases();
		initStatistics();





		if(_distanceType == 0){

			_stats = new SimkaStatistics(_nbBanks, _nbNewBanks, _computeSimpleDistances, _computeComplexDistances, _kmerSize, _sketchSize);
			simka2_loadStatInfos(_database._dir, _database._uniqKmerSpectrumDirs, _database._entries, _kmerSpectrumDirs, _stats, _database._entriesInfos);
			delete _stats;

			DistinctMergedKmerCounter<span> distinctMergedKmerCounter(_partitionId, _kmerSpectrumDirs, _idToOrder);
			distinctMergedKmerCounter.execute();
			u_int64_t nbDistinctMergedKmers = distinctMergedKmerCounter._nbDistinctMergedKmers;
			//string filename = _dirMatrixParts + "/" + Stringify::format("%i", _partitionId) + "-nbDistinctMergedKmers.bin";
			string filename = _dirMatrixParts + "/" + "nbDistinctMergedKmers.bin";
			ofstream file(filename.c_str(), std::ios::binary);
		    file.write((char const*)(&nbDistinctMergedKmers), sizeof(nbDistinctMergedKmers));
		    file.close();

		}
		else if(_distanceType == 1){
			distance();
		}


		//distance();


		//saveStats();

		//delete _stats;
		//delete _processor;


	}

	void parseArgs(){

    	_nbCores =  getInput()->getInt(STR_NB_CORES);
    	_kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	_partitionId = getInput()->getInt(STR_SIMKA2_PARTITION_ID);
    	_nbParitionToProcess = getInput()->getInt(STR_SIMKA2_NB_PARTITION);
    	_databaseDir =  getInput()->getStr(STR_SIMKA2_DATABASE_DIR);
    	_maxDatasets = getInput()->getInt(STR_SIMKA2_DISTANCE_MAX_PROCESSABLE_DATASETS);
    	_nbDistinctMergedKmers = getInput()->getInt(STR_SIMKA2_NB_DISTINCT_MERGED_KMERS);

		_sketchSize = getInput()->getInt(STR_SIMKA_SKETCH_SIZE);
		_distanceType = getInput()->getInt(STR_SIMKA2_DISTANCE_TYPE);

    	_dirMatrixParts = _databaseDir + "/distance/temp_parts";
    	//_dirMatrixMatrixBinary = _databaseDir + "/distance/matrix_binary";

    	//_sketchSize = 10000; //TODO
    	//cout << _partitionId << " " << _nbParitionToProcess << " " << _nbDistinctMergedKmers << endl;
	}

	void createDatabases(){

		_database = Simka2Database(_databaseDir, _maxDatasets);

		_nbBanks = _database._entries.size();
		_nbNewBanks = _nbBanks - _database._nbProcessedDataset;
		//_nbNewBanks = min(_nbNewBanks, _maxDatasets);

		cout << "Database size: " << _nbBanks << endl;
		cout << "Nb new banks to process: " << _nbNewBanks << endl;

		for(size_t i=0; i<_database._entries.size(); i++){
			_idToOrder[_database._entries[i]] = i;
		}

		//cout << "Simka2Distance (createDatabase):   " << _nbBanks << "   " << _nbNewBanks << endl;
		//cout << "uniqDirs: " << _database._uniqKmerSpectrumDirs.size() << endl;
		/*
		_allIds.push_back(datasetID);
		_allKmerSpectrumDirs.push_back(line);

		_newIds.push_back(datasetID);
		_nbBanks = _allIds.size();
		_nbNewBanks = _newIds.size();*/

	}



	void initStatistics(){



		//for(size_t i=0; i<_nbBanks; i++){
		//	cout << _stats-> << endl;
		//}
	}

	void distance(){

		DatasetMergerDistance<span> datasetMergerDistance(_nbBanks, _partitionId, _nbParitionToProcess, _kmerSpectrumDirs, _idToOrder, _nbDistinctMergedKmers, _database, _computeSimpleDistances, _computeComplexDistances, _nbNewBanks, _kmerSize, _sketchSize, _kmerSpectrumDirs, _dirMatrixParts);
		datasetMergerDistance.execute();

		_processor->end();

	}





};











class Simka2Distance: public Tool{
public:

	string _execFilename;

	Simka2Distance(string execFilename): Tool ("Simka2-Distance"){
		_execFilename = execFilename;

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for merged kmer spectrum", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DATABASE_DIR, "dir path to a simka database", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_PARTITION_ID, "processus ID", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_NB_PARTITION, "number of the partition to process", false));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_MAX_PROCESSABLE_DATASETS, "maximum number of datasets that can be processed", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_TYPE, "0: count nb distinct merged kmers, 1: compute distances", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_NB_DISTINCT_MERGED_KMERS, "number of distinct kmers among N datasets", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_MATRIX_DIR, "input filename of k-mer spectrums for which distances has to be computed", true));

	    IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", true));
	    kmerParser->push_back (new OptionOneParam (STR_SIMKA_SKETCH_SIZE, "number of kmers used to compute distances", true));

	    //parser->getParser(STR_NB_CORES)->setVisible(false);
		parser->push_back(kmerParser);

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

    		Simka2DistanceAlgorithm<span>* algo = new Simka2DistanceAlgorithm<span>(p._props);
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
    	Simka2Distance(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



