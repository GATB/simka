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
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;
	struct kxpcomp { bool operator() (kxp l,kxp r) { return (get<0>(r) < get<0>(l)); } } ;

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

	DatasetMergerDistance(size_t nbBanks, size_t partitionId, vector<string>& datasetToMergeDirs, SimkaStatistics* stats, SimkaCountProcessorSimple<span>* processor, map<string, u_int64_t>& idToOrder):
		DiskBasedMergeSort<span>(partitionId, datasetToMergeDirs, idToOrder, true)
    {


		_nbBanks = nbBanks;
		_abundancePerBank.resize(_nbBanks, 0);
		_solidCounter = new SimkaCounterBuilderMerge(_abundancePerBank);
		_isInit = false;
		_nbBankThatHaveKmer = 0;
		_stats = stats;
		_processor = processor;

    }

	void process(Type& kmer, u_int64_t bankId, u_int64_t abundance){
    	//cout << kmer.toString(61) << " " << bankId <<  " "<< abundance << endl;

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


		/*
		cout << kmer.toString(31) << "     ";
		for(size_t i=0; i<counts.size(); i++){
			cout << counts[i] << " ";
		}
		cout << endl;
		cout << "-----" << endl;*/

		_stats->_nbDistinctKmers += 1;

		if(_stats->_computeComplexDistances || _nbBankThatHaveKmer > 1){

			if(_nbBankThatHaveKmer > 1){
				_stats->_nbSharedKmers += 1;
			}

			_processor->process(this->_partitionId, kmer, counts);

		}
	}


	void end(){
		insert(_lastKmer, _abundancePerBank, _nbBankThatHaveKmer);
		delete _solidCounter;
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
		distance();


		saveStats();


		delete _stats;
		delete _processor;

		//writeFinishSignal();

	}

	void parseArgs(){

    	_nbCores =  getInput()->getInt(STR_NB_CORES);
    	_kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	_partitionId = getInput()->getInt(STR_SIMKA2_PARTITION_ID);
    	_databaseDir =  getInput()->getStr(STR_SIMKA2_DATABASE_DIR);

    	_dirMatrixParts = _databaseDir + "/distance/temp_parts";
    	//_dirMatrixMatrixBinary = _databaseDir + "/distance/matrix_binary";

	}

	void createDatabases(){

		_database = Simka2Database(_databaseDir);

		_nbBanks = _database._entries.size();
		_nbNewBanks = _nbBanks - _database._nbProcessedDataset;

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

		_stats = new SimkaStatistics(_nbBanks, _nbNewBanks, _computeSimpleDistances, _computeComplexDistances);
		simka2_loadStatInfos(_databaseDir, _database._uniqKmerSpectrumDirs, _database._entries, _kmerSpectrumDirs, _stats, _database._entriesInfos);

	}

	void distance(){

		_processor = new SimkaCountProcessorSimple<span> (_stats, _nbBanks, _nbNewBanks, _kmerSize, 0);

		DatasetMergerDistance<span> datasetMergerDistance(_nbBanks, _partitionId, _kmerSpectrumDirs, _stats, _processor, _idToOrder);
		datasetMergerDistance.execute();

		_processor->end();


	}

	void saveStats(){
		string filename = _dirMatrixParts + "/" + Stringify::format("%i", _partitionId) + ".gz";
		_stats->save(filename);
	}

	void writeFinishSignal(){
		string finishFilename = _dirMatrixParts + "/" + Stringify::format("%i", _partitionId) + "-success";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
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
	    parser->push_front (new OptionOneParam (STR_SIMKA2_PARTITION_ID, "number of the partition", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_MATRIX_DIR, "input filename of k-mer spectrums for which distances has to be computed", true));

	    IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", true));

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



