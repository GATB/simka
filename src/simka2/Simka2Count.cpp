/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
#include "../core/SimkaUtils.hpp"
#include "Simka2Utils.hpp"
#include "../minikc/MiniKC.hpp"
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include "../utils/SimkaIoUtils.hpp"
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"




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

















/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2ComputeKmerSpectrumAlgorithm : public Algorithm
{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef tuple<Count, StorageItKmerCount<span>*> KmerCount_It;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
	struct kxpcomp { bool operator() (KmerCount_It& l,KmerCount_It& r) { return (get<0>(r).value < get<0>(l).value); } } ;

	u_int64_t _nbReads;


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

	SimkaPartitionWriter<span>* _partitionWriter;

	Simka2ComputeKmerSpectrumAlgorithm(IProperties* options):
		Algorithm("simka", -1, options)
	{

	}

	void execute(){
		_datasetIDbin = 0;

		parseArgs();
		//layoutInputFilename();
		count();
		saveInfos();

		delete _partitionWriter;

		//writeFinishSignal();
	}

	void parseArgs(){

		_options = getInput();

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
			exit(1);
			/*
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
				std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
				exit(1);
			}*/
		}

		//_outputDirTemp = _outputDirTemp;

		if(!System::file().doesExist(_outputDirTemp)){
			std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
			exit(1);
			/*
			int ok = System::file().mkdir(_outputDirTemp, -1);
			if(ok != 0){
				std::cerr << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
				exit(1);
			}*/
		}

		//_outputDirTemp = System::file().getRealPath(_outputDirTemp) + "/";
		//cout << _outputDirTemp << endl;
		//_outputDirTemp += "/" + _datasetID + "_temp" + "/";
		//System::file().mkdir(_outputDirTemp, -1);

		//_options->setStr(STR_URI_OUTPUT_TMP, _outputDirTemp);
		//System::file().mkdir(_outputDirTemp + "/input/", -1);

	}


	/*
	void layoutInputFilename(){

		//string inputDir = _outputDirTemp;
		//ifstream inputFile(_inputFilename.c_str());

		//string bankId = Stringify::format("%i", _datasetID);
		//_banksInputFilename =  inputDir + "__input_simka__"; //_inputFilename + "_dsk_dataset_temp__";
		//IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");
		_banksInputFilename = _outputDirTemp + "/" + _datasetID + "_input";
		//cout << subBankFilename << endl;
		IFile* subBankFile = System::file().newFile(_banksInputFilename, "wb");

		string linePairedDatasets = _inputFilename;
		//cout << linePairedDatasets << endl;
		string linePart;
		vector<string> lineIdDatasets;
		vector<string> linepartPairedDatasets;
		vector<string> linepartDatasets;

		string bankFileContents = "";

		u_int64_t lineIndex = 0;

		stringstream linePairedDatasetsStream(linePairedDatasets);
		while(getline(linePairedDatasetsStream, linePart, ';')){
			linepartPairedDatasets.push_back(linePart);
		}


		//cout << subBankFile->getPath() << endl;
		string subBankContents = "";
		//_nbBankPerDataset.push_back(linepartPairedDatasets.size());
		_nbBankPerDataset = linepartPairedDatasets.size();

		for(size_t i=0; i<linepartPairedDatasets.size(); i++){
			string lineDatasets = linepartPairedDatasets[i];

			linepartDatasets.clear();

			stringstream lineDatasetsStream(lineDatasets);
			while(getline(lineDatasetsStream, linePart, ',')){
				//cout << linePart << endl;
				linepartDatasets.push_back(linePart);
				//cout << "\t" << linePart << endl;
			}

			//bankFileContents += linepartDatasets[0] + "\n";


			for(size_t i=0; i<linepartDatasets.size(); i++){
				string filename = linepartDatasets[i];
				//cout << "lol" << endl;
				//cout << filename << endl;
				//cout << filename.size() << endl;
				//cout << filename << endl;
				subBankContents +=  filename + "\n";
			}

		}

		subBankContents.erase(subBankContents.size()-1);
		subBankFile->fwrite(subBankContents.c_str(), subBankContents.size(), 1);
		subBankFile->flush();
		delete subBankFile;

			//bankFileContents += inputDir + "/" + bankId + "\n";
			//lineIndex += 1;

			//_bankNames.push_back(bankId);


		//}


		//inputFile.close();

		//bankFileContents.erase(bankFileContents.size()-1);
		//bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);
		//bankFile->flush();
		//delete bankFile;

	}*/

	void count(){

		_partitionWriter = new SimkaPartitionWriter<span>(_outputDir, _nbPartitions);

		//vector<string> outInfo;

		IBank* bank = Bank::open(_inputFilename);
		LOCAL(bank);



		u_int64_t nbReads = 0;
		//string tempDir = _outputDirTemp;
		//System::file().mkdir(tempDir, -1);


		//cout << i << endl;
		//string outputDir = p.outputDir + "/comp_part" + to_string(p.datasetId) + "/";

		//cout << "\tinput: " << p.outputDir + "/input/" + p.bankName << endl;

		SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);

		IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _maxNbReads, _nbBankPerDataset);
		// = new SimkaPotaraBankFiltered(bank)
		LOCAL(filteredBank);
		//LOCAL(bank);



		if(_kmerSize <= 15){
			//cerr << "Mini Kc a remettre" << endl;
			//exit(1);
			SimkaMiniKmerCounter<span> miniKc(_options, _kmerSize, filteredBank, _outputDir, _nbPartitions, _abundanceThreshold.first, _abundanceThreshold.second, _partitionWriter);
			miniKc.execute();

			_nbReads = miniKc._nbReads;
		}
		else{

			{

				_h5Filename = _outputDirTemp + "/" +  _datasetID + ".h5";
				Storage* solidStorage = StorageFactory(STORAGE_HDF5).create (_h5Filename, true, false);
				LOCAL(solidStorage);

				ConfigurationAlgorithm<span> configAlgo(filteredBank, _options);
				configAlgo.execute();
				RepartitorAlgorithm<span> repart (filteredBank, solidStorage->getGroup(""), configAlgo.getConfiguration());
				repart.execute ();
				Repartitor* repartitor = new Repartitor();
				LOCAL(repartitor);
				repartitor->load(solidStorage->getGroup(""));


				SimkaAbundanceProcessor<span>* proc = new SimkaAbundanceProcessor<span>(_abundanceThreshold.first, _abundanceThreshold.second);
				CountProcessorDump<span>* procDump = new CountProcessorDump     <span> (solidStorage->getGroup("dsk"), _kmerSize);

				//SimkaCompressedProcessor<span>* proc = new SimkaCompressedProcessor<span>(bags, caches, cacheIndexes, p.abundanceMin, p.abundanceMax);
				std::vector<ICountProcessor<span>* > procs;
				//procs.push_back(proc);
				ICountProcessor<span>* result = 0;

				result = new CountProcessorChain<span> (
						proc,
						procDump,
						NULL
				);


				result->setName ("dsk");


				//SortingCountAlgorithm<span> algo (filteredBank, _options);
				SortingCountAlgorithm<span> algo (filteredBank, configAlgo.getConfiguration(), repartitor,
						procs,
						_options);
				algo.addProcessor(result);
				algo.execute();

				//delete proc;
				//delete procDump;
				//delete result;
				_nbReads = algo.getInfo()->getInt("seq_number");
				_localNbPartitions = algo.getConfig()._nb_partitions;
			}

			partitionKmerCounts();
		}



	}

	void partitionKmerCounts(){



		Storage* solidStorage = 0;
		//string solidsName = _outputDir + "/solid/" +  p.bankName + ".h5";
		//cout << solidsName << endl;
		//bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
		solidStorage = StorageFactory(STORAGE_HDF5).load(_h5Filename);
		LOCAL(solidStorage);

		vector<StorageItKmerCount<span>*> its;
		Partition<Count>& solidKmers = solidStorage->getGroup("dsk").getPartition<Count>("solid");
		for(size_t i=0; i<_localNbPartitions; i++){
			//cout << "lala" << endl;
			//cout << solidKmers[i].iterable()->getNbItems () << endl;
			its.push_back(new StorageItKmerCount<span>(solidKmers[i].iterable()->iterator()));
		}

		std::priority_queue< KmerCount_It, vector<KmerCount_It>,kxpcomp > pq;
		//StorageItKmerCount<span>* bestIt;
		size_t bestPart;

		for(size_t i=0; i<its.size(); i++){
			its[i]->_it->first();
			//StorageItKmerCount<span>* it = its[i];
			//it->_it->first();
		}

		//fill the  priority queue with the first elems
		for (size_t ii=0; ii<its.size(); ii++)
		{
			//pq.push(Kmer_BankId_Count(ii,its[ii]->value()));
			if (!its[ii]->_it->isDone()){
				pq.push(KmerCount_It(its[ii]->item(), its[ii]));
			}
		}

		StorageItKmerCount<span>* bestIt;

		if (pq.size() != 0) // everything empty, no kmer at all
		{
			//get first pointer
			//bestPart =
			//bestIt = get<3>(pq.top()); pq.pop();
			KmerCount_It kmerCountIt = pq.top(); pq.pop();
			_partitionWriter->insert(get<0>(kmerCountIt).value, _datasetIDbin, get<0>(kmerCountIt).abundance);

			bestIt = get<1>(kmerCountIt);


			while(1){

				if (! bestIt->next())
				{
					//reaches end of one array
					if(pq.size() == 0){
						break;
					}

					//otherwise get new best
					//best_p = get<1>(pq.top()) ; pq.pop();
					bestIt = get<1>(pq.top()); pq.pop();
				}

				pq.push(KmerCount_It(bestIt->item(), bestIt));
				//pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

				bestIt = get<1>(pq.top()); pq.pop();

				_partitionWriter->insert(bestIt->item().value, _datasetIDbin, bestIt->item().abundance);
				/*
				//_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
				size_t part = oahash(bestIt->item().value) % _nbPartitions;
				//cout << part << " " << p.nbPartitions << endl;
				//cout << bestIt->item().value.toString(21) << " " << bestIt->item().abundance << endl;
				cachedBags[part]->insert(Kmer_BankId_Count(bestIt->item().value, _datasetIDbin, bestIt->item().abundance));

				//Kmer_BankId_Count item(kmer, _bankIndex, count[0]);
				_nbDistinctKmerPerParts[part] += 1;
				_nbKmerPerParts[part] += bestIt->item().abundance;
				_chordNiPerParts[part] += pow(bestIt->item().abundance, 2);*/
			}
		}

		_partitionWriter->end();

		for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}
		//System::file().remove();
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

        ofstream outputInfoFile(_outputDir + "/merge_info.bin", std::ios::binary);

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

	void writeFinishSignal(){

		string finishFilename = _outputDir + "/success";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
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

    		Simka2ComputeKmerSpectrumAlgorithm<span>* algo = new Simka2ComputeKmerSpectrumAlgorithm<span>(p._props);
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



