

#include <gatb/gatb_core.hpp>
#include "ProbabilisticDictionary.hpp"
#include "SimkaMinUtils.hpp"

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


typedef ProbabilisticDictionary<u_int16_t, u_int16_t> ProbabilisticDict;


































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
	ProbabilisticDict* _selectedKmersIndex;

	bool _isMaster;
	bool _exists;

	CountKmerCommand(size_t kmerSize, ProbabilisticDict* kmerIndex)
	: _model(kmerSize), _itKmer(_model)
	{
		_selectedKmersIndex = kmerIndex;
		_kmerSize = kmerSize;
		_isMaster = true;
	}

	CountKmerCommand(const CountKmerCommand& copy)
	: _model(copy._kmerSize), _itKmer(_model)
	{
		_kmerSize = copy._kmerSize;
		_isMaster = false;
		_selectedKmersIndex = copy._selectedKmersIndex;
	}


	void operator()(Sequence& sequence){
		//cout << sequence.getIndex() << endl;
		//cout << sequence.toString() << endl;

		_itKmer.setData(sequence.getData());

		for(_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			//cout << _itKmer->value().toString(_kmerSize) << endl;

			u_int64_t kmerValue = _itKmer->value().getVal();
			u_int64_t index = _selectedKmersIndex->getIndex(kmerValue, _exists);
			//if(index > 1000){
			//	cout << kmerValue << " " << index  << endl;
			//}
			if(_exists){
				_selectedKmersIndex->increment(index);
				//cout << kmerValue << " " << index  << endl;
			}

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
class SimkaMinAlgorithm : public Algorithm
{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type  KmerType;
	typedef typename Kmer<span>::ModelCanonical ModelCanonical;
	typedef typename Kmer<span>::ModelCanonical::Iterator ModelCanonicalIterator;


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
	//vector<u_int64_t> _nbReadsPerDataset;

	//string _outputFilenameSuffix;

	//u_int64_t _totalKmers;
	vector<size_t> _nbBankPerDataset;
	u_int64_t _nbBitPerKmers;
	u_int64_t _bloomFilterMaxMemoryBits;

	//string _largerBankId;
	//bool _computeSimpleDistances;
	//bool _computeComplexDistances;
	//bool _keepTmpFiles;


	ProbabilisticDict* _selectedKmersIndex;
	Bloom<Type>*  _bloomFilter;


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

		computeNbUsedKmers();
		selectKmers();

		//reprise: la stratégie actuelle est trop lente car pas assez e calcul par rapport au decompression gz. abandonner l'idée du dispatcher et paralleliser au niveau des jeux de données
		string filename = _outputDirTemp + "/selectedKmers.bin";
		file_binary<u_int64_t> keyFile(filename.c_str());
		_selectedKmersIndex = new ProbabilisticDict(_nbUsedKmers, keyFile, _nbCores);

		countKmers();
	}

	void parseArgs(){
		//_keepTmpFiles = _options->get(STR_SIMKA_KEEP_TMP_FILES);
		_nbBitPerKmers = 12;
		_maxMemory = _args->getInt(STR_MAX_MEMORY);
		_nbCores = _args->getInt(STR_NB_CORES);
		_inputFilename = _args->getStr(STR_URI_INPUT);
		_outputDir = _args->get(STR_URI_OUTPUT) ? _args->getStr(STR_URI_OUTPUT) : "./";
		_outputDirTemp = _args->get(STR_URI_OUTPUT_TMP) ? _args->getStr(STR_URI_OUTPUT_TMP) : "./";
		_kmerSize = _args->getInt(STR_KMER_SIZE);
		_nbUsedKmers = _args->getInt(STR_SIMKA_NB_KMERS_USED);
		//_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
		//_abundanceThreshold.second = min((u_int64_t)_options->getInt(STR_KMER_ABUNDANCE_MAX), (u_int64_t)(999999999));

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


	void computeNbUsedKmers(){
		cout << endl << endl;
		u_int64_t bloomFilterMaxMemoryByte = _maxMemory*GBYTE;
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
		_bloomFilter = new BloomCacheCoherent<Type> (_bloomFilterMaxMemoryBits, 7);

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

	    for (size_t i=0; i<itBanks.size(); i++)
	    {

	    	Iterator<Sequence>* itSeq = itBanks[i];
	    	u_int64_t nbKmerInserted = 0;
	    	bool isDone = false;

			for(itSeq->first(); !itSeq->isDone(); itSeq->next()){
				Sequence& sequence = itSeq->item();

				itKmer.setData(sequence.getData());

				for(itKmer.first(); !itKmer.isDone(); itKmer.next()){

					Type kmer = itKmer->value();
					if(kmer.getVal() == 0) continue; //todo can be optimized maybe ?

					if(_bloomFilter->contains(kmer)){

					}
					else{
						u_int64_t kmerValue = kmer.getVal();
						selectKmersFile.write((const char*)&kmerValue, sizeof(kmerValue));
						_bloomFilter->insert(kmer);
						nbKmerInserted += 1;
						nbKmerInsertedTotal += 1;
						if(nbKmerInserted >= _nbUsedKmersPerDataset){
							isDone = true;
							break;
						}
					}
				}

				if(isDone){
					cout << nbKmerInserted << endl;
					break;
				}

			}

			itBanks[i]->finalize();
	    }

	    _nbUsedKmers = nbKmerInsertedTotal;
	    cout << "Final number of selected kmers: " << _nbUsedKmers << endl;
	    selectKmersFile.close();
	    delete _bloomFilter;
	}

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

	void countKmers(){



		IBank* bank = Bank::open(_banksInputFilename);
		LOCAL(bank);
		Iterator<Sequence>* it = bank->iterator();
		std::vector<Iterator<Sequence>*> itBanks =  it->getComposition();

		size_t threadId = 0;
		_maxRunningThreads = _nbCores;
		//vector<thread> threads; //(_nbCores);
		//_isThreadRunning = vector<bool>(_nbCores);
		_nbRunningThreads = 0;

		for (size_t i=0; i<itBanks.size(); i++){
			//cout << "processing bank: " << i << endl;

			Iterator<Sequence>* itSeqSimka = new SimkaInputIterator<Sequence> (itBanks[i], _nbBankPerDataset[i], _maxNbReads);
			//LOCAL(itSeqSimka);

			//size_t threadId = obtainThreadId();

			thread* t = new thread(&SimkaMinAlgorithm::countKmersOfDataset, this,  itSeqSimka, threadId);
			_threads.push_back(t);
			_runningThreadIds.push_back(threadId);
			threadId += 1;
			_nbRunningThreads += 1;
			//_isThreadRunning[threadId] = true;
			//_nbRunningThreads[i] += 1;

			if(_nbRunningThreads >= _maxRunningThreads){
				countKmersWait();
			}



			//---------------------------
			//todo !!!!!!!!!!!!!!!!!!! delete itSeqSimka, thread, itBanks[i]->finalize();
			//---------------------------
		}

		countKmersJoinThreads();

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

	void countKmersOfDataset(Iterator<Sequence>* itSeqSimka, size_t threadId){

		countKmersMutex.lock();
		cout << "starting thread: " << threadId << endl;
		size_t nbCores = 1;
		IDispatcher* dispatcher = new Dispatcher (nbCores);
		countKmersMutex.unlock();

		dispatcher->iterate(itSeqSimka, CountKmerCommand<span>(_kmerSize, _selectedKmersIndex), 1000);

		countKmersMutex.lock();
		cout << "\t thread " <<  threadId << " End" << endl;
		_finishedThreads.push_back(threadId);
		countKmersMutex.unlock();
	}

	void countKmersWait(){
        while(1){


			countKmersMutex.lock();

			bool isThreadAvailbale = false;


			for(size_t i=0; i<_finishedThreads.size(); i++){
				size_t threadId = _finishedThreads[i];

				//_runningThreadIds.erase(std::remove(_runningThreadIds.begin(), _runningThreadIds.end(), threadId), _runningThreadIds.end());
				auto it = find(_runningThreadIds.begin(), _runningThreadIds.end(), threadId);
				int pos = distance(_runningThreadIds.begin(), it);

				cout << "\t removing thread " <<  threadId << " (pos: "  << pos << ")" << endl;

				_runningThreadIds.erase(_runningThreadIds.begin()+pos);
				_threads[pos]->join();
				delete _threads[pos];
				_threads.erase(_threads.begin()+pos);

				_nbRunningThreads -= 1;
				isThreadAvailbale = true;

			}


			if(isThreadAvailbale){
				_finishedThreads.clear();
				cout << _runningThreadIds.size() << " " << _threads.size() << endl;
				countKmersMutex.unlock();
				break;
			}

			countKmersMutex.unlock();

			sleep(1);

        }

	}

	void countKmersJoinThreads(){
        while(_nbRunningThreads > 0)
        	countKmersWait();
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
	    readParser->push_back (new OptionOneParam (STR_SIMKA_MAX_READS.c_str(), "maximum number of reads per sample to process. Can be -1: use all reads. Can be 0: estimate it", false, "-1" ));
	    //readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE.c_str(), "minimal size a read should have to be kept", false, "0" ));
	    //readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(), "minimal Shannon index a read should have to be kept. Float in [0,2]", false, "0" ));

	    //Core parser
	    IOptionsParser* coreParser = new OptionsParser ("core");
	    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
	    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (GB)", false, "8"));
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



