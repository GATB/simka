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

#include "SimkaAlgorithm.hpp"

static const char* strProgressPartitionning = "Simka: Step 1: partitioning    ";
static const char* strProgressCounting =      "Simka: Step 2: counting kmers  ";


template<size_t span>
SimkaCountProcessor<span>::SimkaCountProcessor (SimkaStatistics& stats, size_t nbBanks, const pair<size_t, size_t>& abundanceThreshold, SIMKA_SOLID_KIND solidKind, bool soliditySingle, IteratorListener* progress, const SimkaDistanceParam& distanceParams) :
_progress(progress), _stats(stats)
{

	// We configure the vector for the N.(N+1)/2 possible pairs
	//_countTotal.resize (_nbBanks*(_nbBanks+1)/2);

	_nbBanks = nbBanks;
	_abundanceThreshold = abundanceThreshold;
	_solidKind = solidKind;
	_soliditySingle = soliditySingle;

	_localStats = new SimkaStatistics(_nbBanks);
	_distanceParams = distanceParams;

	_nbKmerCounted = 0;
}

template<size_t span>
SimkaCountProcessor<span>::~SimkaCountProcessor () {

	if(_progress){ //Simka_min
		if(_nbKmerCounted > 0){
			_progress->inc(_nbKmerCounted);
			_nbKmerCounted = 0;
		}
	}

	delete _localStats;
}


template<size_t span>
void SimkaCountProcessor<span>::finishClones (std::vector<ICountProcessor<span>*>& clones)
{
	for (size_t i=0; i<clones.size(); i++){

		if (SimkaCountProcessor* clone = dynamic_cast<SimkaCountProcessor*> (clones[i])){
			finishClone(clone);
			//for (size_t i=0; i<this->_countTotal.size(); i++)  { this->_countTotal[i] += clone->_countTotal[i];  }
		}
	}
}

template<size_t span>
void SimkaCountProcessor<span>::finishClone(SimkaCountProcessor<span>* clone){
	//cout << _stats << "   " << &clone->_stats << endl;
	_stats += *clone->_localStats;
}

template<size_t span>
bool SimkaCountProcessor<span>::isSolidVector(const CountVector& counts){

	//if(_solidKind == SIMKA_SOLID_KIND::RANGE){
	//}

	//bool isSolid_ = false;

	for(size_t i=0; i<counts.size(); i++){

		//cout << counts[i] << " " << _abundanceThreshold.first << endl;
		if(counts[i] >= _abundanceThreshold.first)
			return true;

		/*
		CountNumber abundance = counts[i];

		if(abundance == 0) continue;

		if(isSolid(abundance)){
			return true;
		}*/

		//nbBanks += 1;
		//if(nbBanks > 1){
		//	isSolid = true;
		//	break;
		//}
	}

	return false;

}


template<size_t span>
bool SimkaCountProcessor<span>::isSolid(CountNumber count){
	return count >= _abundanceThreshold.first && count <= _abundanceThreshold.second;
}


template<size_t span>
bool SimkaCountProcessor<span>::process (size_t partId, const Type& kmer, const CountVector& counts, CountNumber sum){


	if(_progress){ //Simka_min
		if(_nbKmerCounted > 500000){
			_progress->inc(_nbKmerCounted);
			_nbKmerCounted = 0;
		}
	}

	//return false;

	_totalAbundance = 0;
	_localStats->_nbDistinctKmers += 1;

	for(size_t i=0; i<counts.size(); i++){

		CountNumber abundance = counts[i];
		//_nbKmerCounted += abundance;
		//_stats._speciesAbundancePerDataset[i].push_back(abundance);

		//cout << counts[i] << " ";
		_localStats->_nbKmers += abundance;
		_localStats->_nbKmersPerBank[i] += abundance;
		_totalAbundance += abundance;
	}

	//cout << kmer.toString(31) << endl;

	//cout << endl;

	if(_progress){ //Simka_min
		_localStats->_nbSolidKmers += 1;
		computeStats(counts);
	}
	else{

		if(_abundanceThreshold.first > 1){
			if(!isSolidVector(counts))
				return false;
		}

		_localStats->_nbSolidKmers += 1;

		if(_soliditySingle){
			CountVector counts2(counts);
			for(size_t i=0; i<counts.size(); i++){
				//if(!isSolid(counts[i]))
				if(counts[i] < _abundanceThreshold.first)
					counts2[i] = 0;
			}
			computeStats(counts2);
		}
		else{
			computeStats(counts);
		}

	}

	/*
	CountVector counts2;
	for(int i=0; i<counts.size(); i++){
		if(counts[i] < _abundanceThreshold.first)
			counts2.push_back(0);
		else
			counts2.push_back(counts[i]);
	}*/


	return true;
}


template<size_t span>
void SimkaCountProcessor<span>::computeStats(const CountVector& counts){

	int nbBanksThatHaveKmer = 0;
	//u_int64_t totalAbundance = 0;



	for(size_t i=0; i<counts.size(); i++){

		CountNumber abundanceI = counts[i];

		if(abundanceI){
			nbBanksThatHaveKmer += 1;
			_localStats->_nbSolidDistinctKmersPerBank[i] += 1;
			_localStats->_nbSolidKmersPerBank[i] += abundanceI;
		}

		u_int64_t bc;

		for(size_t j=i+1; j<counts.size(); j++){
			CountNumber abundanceJ = counts[j];

			/*
			if(abundanceI < abundanceJ){
				_localStats->_brayCurtisNumerator[i][j] += abundanceI;
				_localStats->_brayCurtisNumerator[j][i] += abundanceI;
			}
			else{
				_localStats->_brayCurtisNumerator[i][j] += abundanceJ;
				_localStats->_brayCurtisNumerator[j][i] += abundanceJ;
			}*/
			//updateBrayCurtis(i, abundanceI, j, abundanceJ);

			if(_distanceParams._computeBrayCurtis){
				bc = min(abundanceI, abundanceJ);
				_localStats->_brayCurtisNumerator[i][j] += bc;
				_localStats->_brayCurtisNumerator[j][i] += bc;
			}

			if(abundanceI && abundanceJ){
				_localStats->_matrixNbSharedKmers[i][j] += abundanceI;
				_localStats->_matrixNbSharedKmers[j][i] += abundanceJ;
				_localStats->_matrixNbDistinctSharedKmers[i][j] += 1;
				_localStats->_matrixNbDistinctSharedKmers[j][i] += 1;
				//updateKullbackLeibler(i, abundanceI, j, abundanceJ);
			}

		}

	}

	_localStats->_nbDistinctKmersSharedByBanksThreshold[nbBanksThatHaveKmer-1] += 1;
	_localStats->_nbKmersSharedByBanksThreshold[nbBanksThatHaveKmer-1] += _totalAbundance;

	if(_totalAbundance == 1){
		//if( == 1){
		_localStats->_nbErroneousKmers += 1;
		//}
	}
	//else if(nbBanksThatHaveKmer == counter.size()){
	//}


}


template<size_t span>
void SimkaCountProcessor<span>::updateBrayCurtis(int bank1, CountNumber abundance1, int bank2, CountNumber abundance2){
	/*
    //n = len(X)
    bc_num = 0
    bc_den = 0
    for i in range(n):
        if (X[i] + Y[i] > 0):
            bc_num += abs(abundance1-abundance2)
            bc_den += abundance1 + abundance2
    bc = bc_num/bc_den
    #
    return bc*/


	//_localStats->_brayCurtisNumerator[bank1][bank2] += abs(abundance1-abundance2);
	_localStats->_brayCurtisNumerator[bank1][bank2] += min(abundance1, abundance2);
	_localStats->_brayCurtisNumerator[bank2][bank1] += min(abundance1, abundance2);
}
































template<size_t span>
SimkaAlgorithm<span>::SimkaAlgorithm(IProperties* options)
:
Algorithm("simka", -1, options),
_progress (0), _tmpPartitionsStorage(0), _tmpPartitions(0)
{


	_stats = 0;
	//_simkaDistance = 0;
	_banks = 0;
	_processor = 0;

	_options = options;


	_maxMemory = getInput()->getInt(STR_MAX_MEMORY);
    _nbCores = getInput()->getInt(STR_NB_CORES);
	_inputFilename = _options->getStr(STR_URI_INPUT);
	_outputDir = _options->get(STR_URI_OUTPUT) ? _options->getStr(STR_URI_OUTPUT) : "./";
	_outputDirTemp = _options->get(STR_URI_OUTPUT_TMP) ? _options->getStr(STR_URI_OUTPUT_TMP) : "./";
	_kmerSize = _options->getInt(STR_KMER_SIZE);
	_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
	_abundanceThreshold.second = _options->getInt(STR_KMER_ABUNDANCE_MAX);
	_soliditySingle = _options->get(STR_SIMKA_SOLIDITY_PER_DATASET);
	_nbMinimizers = getInput()->getInt(STR_KMER_PER_READ);
	//_maxDisk = getInput()->getInt(STR_MAX_DISK);

	//read filter
	_maxNbReads = _options->getInt(STR_SIMKA_MAX_READS);
	_minReadSize = _options->getInt(STR_SIMKA_MIN_READ_SIZE);
	_minReadShannonIndex = _options->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
	_minReadShannonIndex = std::max(_minReadShannonIndex, 0.0);
	_minReadShannonIndex = std::min(_minReadShannonIndex, 2.0);


	//string maxDisk = "";
	//if(_options->get(STR_MAX_DISK)){
	//	maxDisk = _options->getStr(STR_MAX_DISK);
	//	cout << maxDisk << endl;
	//}
	//_multiStorage = new MultiDiskStorage<Type>(_options->getStr(STR_URI_OUTPUT_DIR), _options->getStr(STR_MAX_DISK));

	   // vector<string> _tempDirMaxDisk

	_totalKmers = 0;
	/*
	if(_options->getInt(STR_VERBOSE) != 0){
		cout << "Filter options" << endl;
		cout << "\tMax reads per dataset:  " <<  _maxNbReads << endl;
		cout << "\tMin read size:  " <<  _minReadSize << endl;
		cout << "\tMin Shannon index:  " <<  _minShannonIndex << endl;
	}*/
	//if(_maxNbReads == 0)
	//	_maxNbReads = -1;

	//cout << _maxNbReads << endl;
	//cout << _soliditySingle << endl;
	/*
	string solidKindStr = _options->getStr(STR_SOLIDITY_KIND);
	if(solidKindStr == "range"){
		_solidKind = SIMKA_SOLID_KIND::RANGE;
	}
	else if(solidKindStr == "sum"){
		_solidKind = SIMKA_SOLID_KIND::SUM;
	}

	cout << solidKindStr << " " << solidKindStr << endl;*/
	//_kmerSize = _options->get(STR_KMER_SIZE) ? _options->getInt(STR_KMER_SIZE) : 31;
	//_abundanceMin = _options->get(STR_KMER_ABUNDANCE_MIN) ? _options->getInt(STR_KMER_ABUNDANCE_MIN) : 0;
	//_maxMemory = props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 2000;
	//_outputTempDir = props->get(STR_URI_OUTPUT_DIR) ? props->getStr(STR_URI_OUTPUT_DIR) : System::file().getDirectory(_inputFilename);
	//_outputFilename = props->get(STR_URI_OUTPUT) ? props->getStr(STR_URI_OUTPUT) : System::file().getDirectory(_inputFilename) + "/" + System::file().getBaseName(_inputFilename) + "_output.fasta";
	//_nbCores = getInput()->getInt(STR_NB_CORES);

	//cout << "Input filename: " << _inputFilename << endl;
	//cout << "Kmer size: " << _kmerSize << endl;
	//cout << "Abundance min: " << _abundanceMin << endl;
	//cout << "Max memory: " << _maxMemory << endl;
	//cout << "Output temp dir: " << _outputTempDir << endl;
	//cout << "Output filename: " << _outputFilename << endl;


	_banksInputFilename = _inputFilename + "_dsk_dataset_temp__";


}

template<size_t span>
SimkaAlgorithm<span>::~SimkaAlgorithm() {
}



template<size_t span>
void SimkaAlgorithm<span>::execute() {

	if(!System::file().doesExist(_outputDir)){
		int ok = System::file().mkdir(_outputDir, -1);
		if(ok != 0){
	        std::cout << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
	        return;
		}
	}

	if(_nbMinimizers > 0){
		executeSimkamin();
		return;
	}

	layoutInputFilename();
	createBank();

	/*
#ifdef SIMKA_FUSION

	cout << "SimkaAlgo.cpp 1" << endl;
	clear();
	delete _banks;
	cout << "SimkaAlgo.cpp 2" << endl;
	SimkaFusion<span>* simkaFusion = new SimkaFusion<span>(_options, _inputFilename, _outputDir, _outputDirTemp, _nbReadsPerDataset, _maxNbReads);
	simkaFusion->execute();
	return;
#endif*/

	count();

	outputMatrix();
	//outputHeatmap();

	if(_options->getInt(STR_VERBOSE) != 0){
		_stats->print();
		print();
	}

	clear();
}

template<size_t span>
std::vector<size_t> SimkaAlgorithm<span>::getNbCoresList()
{
    std::vector<size_t> result;

    for (size_t p=0; p<_nbPartitions; )
    {
        u_int64_t ram_total = 0;
        size_t i=0;
        for (i=0; i< _nbCores && p<_nbPartitions
            && (ram_total ==0  || ((ram_total+(_nbKmerPerPartitions[p]*getSizeofPerItem()))  <= _maxMemory*MBYTE)) ; i++, p++)
        {
            ram_total += _nbKmerPerPartitions[p]*getSizeofPerItem();
        }

        result.push_back (i);
    }

    return result;
}


template<size_t span>
void SimkaAlgorithm<span>::executeSimkamin() {

	layoutInputFilename();
	createBank();











	//_stats = new SimkaStatistics(_nbBanks);
	//SortingCountAlgorithm<span> sortingCount (_banks, _options);
	//_processor = new SimkaCountProcessor<span> (*_stats, _nbBanks, _abundanceThreshold, _solidKind, _soliditySingle);
	//_processor->use();
	u_int64_t nbReadToProcess = _maxNbReads * _nbBanks;
	u_int64_t maxNbKmers = nbReadToProcess * _nbMinimizers;
	u_int64_t maxNbKmersMemoryB = maxNbKmers * getSizeofPerItem();
	size_t memoryB = _maxMemory * MBYTE;
	size_t perfectOpenFiles = 0;
	size_t max_open_files = System::file().getMaxFilesNumber() / 1.1;

	//cout << maxNbKmers << endl;
	if(maxNbKmersMemoryB < memoryB){
		perfectOpenFiles = _nbCores;
	}
	else{
		size_t memoryPerCoreB = memoryB / _nbCores;
		memoryPerCoreB /= 1.1;
		perfectOpenFiles = maxNbKmersMemoryB / memoryPerCoreB;
	}

	if(perfectOpenFiles > max_open_files){
		perfectOpenFiles = max_open_files;
	}
	perfectOpenFiles = min(perfectOpenFiles, (size_t)2000);
	perfectOpenFiles = max(perfectOpenFiles, _nbCores);

	//cout << perfectOpenFiles << endl;
	//cout << max_open_files << endl;
	_nbPartitions = perfectOpenFiles;
	cout << "Nb partitions: " << _nbPartitions << endl;


	_nbKmerPerPartitions.resize(_nbPartitions, 0);
	_nbk_per_radix_per_part.resize(256);
    for(size_t ii=0; ii<256; ii++)
	{
		_nbk_per_radix_per_part[ii].resize(_nbPartitions, 0);
	}

    //_multiStorage->createRepartition(_nbPartitions, maxNbKmersMemoryB/_nbPartitions);
    //_multiStorage->createStorages();

    string tmpStorageName = _outputDirTemp + "/" + System::file().getTemporaryFilename("dsk_partitions");
    setPartitionsStorage (StorageFactory(STORAGE_FILE).create (tmpStorageName, true, false));
    setPartitions        (0); // close the partitions first, otherwise new files are opened before  closing parti from previous pass
    setPartitions        ( & (*_tmpPartitionsStorage)().getPartition<Type> ("parts", _nbPartitions));
    cout << "Tmp storage: " << tmpStorageName << endl;



    vector<Iterator<Sequence>*> itBanks =  _banks->iterator()->getComposition();

    cout << endl;
    //int total = 0;
	//u_int64_t progressUpdateStep = nbReadToProcess / 100;
	//progressUpdateStep = max(progressUpdateStep, (u_int64_t)10000);
	//cout << progressUpdateStep << endl;
    //size_t nbIterations = (1 + _processors.size()) * _config._volume * MBYTE / sizeof(Type);
    setProgress (new ProgressSynchro (
        createIteratorListener (nbReadToProcess, strProgressPartitionning),
        System::thread().newSynchronizer())
    );
    _progress->init ();

    for (size_t i=0; i<itBanks.size(); i++){

    	//cout << i << endl;

        getDispatcher()->iterate (itBanks[i], FillPartitions<span> (_progress, _nbMinimizers, _nbPartitions, _kmerSize, _maxNbReads, _tmpPartitions, _nbk_per_radix_per_part, _minKmerShannonIndex), 1000, true);

//#ifdef MULTI_DISK
        //        _multiStorage->flush();
        //#else
        _tmpPartitions->flush();
        //#endif

        vector<size_t> nbItems;
        for (size_t p=0; p<_nbPartitions; p++)
        {
        	//#ifdef MULTI_DISK
        	//        u_int64_t nbItem = _multiStorage->getPartition(p).getNbItems();
        	//#else
        u_int64_t nbItem = (*_tmpPartitions)[p].getNbItems();
        //#endif

            nbItems.push_back (nbItem);
            //_nbKmerPerPartitions[p] += nbItem;
            //total += nbItem;
        }
        _nbKmersPerPartitionPerBank.push_back (nbItems);

		//GR: close the input bank here with call to finalize
		itBanks[i]->forget(); //->finalize();
    }


    _totalKmers = 0;
    for (size_t p=0; p<_nbPartitions; p++){
    	//#ifdef MULTI_DISK
    	//    	u_int64_t nbItem = _multiStorage->getPartition(p).getNbItems();
    	//#else
        u_int64_t nbItem = (*_tmpPartitions)[p].getNbItems();
        //#endif

    	_nbKmerPerPartitions[p] = nbItem;
        _totalKmers += nbItem;
    }
    cout << endl << "Total kmers: "  << _totalKmers << endl;



    //_nbCores = 1;








    vector<size_t> coreList = getNbCoresList();

    cout << endl << "Nb cores list:  ";
    for(size_t i=0; i<coreList.size(); i++){
    	cout << coreList[i] << " ";
    }
    cout << endl;


    setProgress (new ProgressSynchro (
        createIteratorListener (_totalKmers, strProgressCounting),
        System::thread().newSynchronizer())
    );
    _progress->init ();


	_stats = new SimkaStatistics(_nbBanks);
	SimkaDistanceParam distanceParams(_options);
    _processor = new SimkaCountProcessor<span> (*_stats, _nbBanks, _abundanceThreshold, _solidKind, _soliditySingle, _progress, distanceParams);
    _processor->use();

    MemAllocator pool (_nbCores);
	u_int64_t memoryPoolSize = _maxMemory*MBYTE;

    size_t p = 0;
    for (size_t i=0; i<coreList.size(); i++)
    {
        vector<ICommand*> cmds;

        vector<ICountProcessor<span>*> clones;

        size_t currentNbCores = coreList[i];

        u_int64_t mem = (_maxMemory*MBYTE)/currentNbCores;

        size_t cacheSize = min ((u_int64_t)(200*1000), mem/(50*sizeof(Count)));

        /** We build a list of 'currentNbCores' commands to be dispatched each one in one thread. */
        for (size_t j=0; j<currentNbCores; j++, p++)
        {
            ISynchronizer* synchro = System::thread().newSynchronizer();
            LOCAL (synchro);

            ICountProcessor<span>* processorClone = _processor->clone ();

            processorClone->use();
            clones.push_back (processorClone);


            uint64_t memoryPartition = (_nbKmerPerPartitions[p]*getSizeofPerItem()); //in bytes
            //cout << _nbKmerPerPartitions[p] << endl;
            //cout << "\t" << getSizeofPerItem() << endl;
            //cout << "\t" << memoryPartition << endl;

            ICommand* cmd = 0;


			if (memoryPartition >= memoryPoolSize)
			{
				static const int EXCEED_FACTOR = 2;

				if (memoryPartition  < EXCEED_FACTOR*memoryPoolSize)
				{
					memoryPoolSize = memoryPartition;
				}
				else
				{
					unsigned long system_mem = System::info().getMemoryPhysicalTotal();
					memoryPoolSize = memoryPartition;

					if (memoryPoolSize > system_mem*0.95)
					{
						throw Exception ("memory issue: %lld bytes required, %lld bytes set by command-line limit, %lld bytes in system memory",
							memoryPartition, memoryPoolSize, system_mem
						);
					}
					else
						cout << "Warning: forced to allocate extra memory: " << memoryPoolSize / MBYTE << " MB" << endl;

				}
			}

		   //if capa pool ==0, reserve max memo , pass pool to partibyvec, will be used  for vec kmers
			if (pool.getCapacity() == 0)  {  pool.reserve (memoryPoolSize); }
			else if (memoryPoolSize > pool.getCapacity()) { pool.reserve(0); pool.reserve (memoryPoolSize); }


			vector<size_t> nbItemsPerBankPerPart;
			for (size_t i=0; i<_nbKmersPerPartitionPerBank.size(); i++)
			{
				nbItemsPerBankPerPart.push_back (_nbKmersPerPartitionPerBank[i][p] - (i==0 ? 0 : _nbKmersPerPartitionPerBank[i-1][p]) );
			}

			//cmd = new PartitionsByVectorCommand<span> (
			//	(*_tmpPartitions)[p], processorClone, cacheSize, _progress, _fillTimeInfo,
			//	pInfo, pass, p, _config._nbCores_per_partition, _config._kmerSize, pool, nbItemsPerBankPerPart
			//);
			int nbCore = 1;
			cmd = new PartitionCommand<span> (
					//_multiStorage->getPartition(p), processorClone, cacheSize, 0, p, nbCore, _kmerSize, pool, nbItemsPerBankPerPart, _nbKmerPerPartitions, _nbk_per_radix_per_part
					(*_tmpPartitions)[p], processorClone, cacheSize, 0, p, nbCore, _kmerSize, pool, nbItemsPerBankPerPart, _nbKmerPerPartitions, _nbk_per_radix_per_part
			);

            cmds.push_back (cmd);

        }

        getDispatcher()->dispatchCommands (cmds, 0);

        _processor->finishClones (clones);
        for (size_t i=0; i<clones.size(); i++)  { delete clones[i]; }  clones.clear();

        //cout << pool.getCapacity()/(double)MBYTE << " " << pool.getUsedSpace()/(double)MBYTE << " " << (pool.getUsedSpace()*100) / (double)pool.getCapacity() << endl;

        pool.free_all();
    }




    //_multiStorage->remove();
    _tmpPartitions->remove();
    //return;
	outputMatrix();
	//outputHeatmap();

	if(_options->getInt(STR_VERBOSE) != 0){
		_stats->print();
		print();
	}

	clear();
}

template<size_t span>
void SimkaAlgorithm<span>::layoutInputFilename(){

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

			bankFileContents += System::file().getBaseName(subBankFilename) + "\n";
			_nbBankPerDataset.push_back(linePartList.size() - 1); //linePartList.size() - 1 = nb sub banks
			//_nbReadsPerDataset.push_back(ceil(_maxNbReads / (float)()));
		}

		lineIndex += 1;
	}

	bankFileContents.erase(bankFileContents.size()-1);
	//bankFileContents.pop_back(); // "remove last /n

	bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);

	inputFile.close();
	//delete inputFile;
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


template<size_t span>
void SimkaAlgorithm<span>::createBank(){

	IBank* bank = Bank::open(_banksInputFilename);
	_nbBanks = bank->getCompositionNb();

	if(_maxNbReads == 0){
		if(_options->getInt(STR_VERBOSE) != 0)
			cout << "-maxNbReads is not defined. Simka will estimating it..." << endl;
		_maxNbReads = bank->estimateNbItems() / _nbBanks;
		_maxNbReads -= (_maxNbReads/10);
		if(_options->getInt(STR_VERBOSE) != 0)
			cout << "Max nb reads: " << _maxNbReads << endl << endl;
	}

	for(size_t i=0; i<_nbBankPerDataset.size(); i++){
		//cout << _maxNbReads << " " << _nbBankPerDataset[i] << endl;
		_nbReadsPerDataset.push_back( ceil(_maxNbReads / (float)(_nbBankPerDataset[i])) );
	}


	SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);


	_banks = new SimkaBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _nbReadsPerDataset, _maxNbReads*_nbBanks);



	//cout << bank->estimateNbItems() << endl;



}

template<size_t span>
void SimkaAlgorithm<span>::count(){

	_stats = new SimkaStatistics(_nbBanks);

	SortingCountAlgorithm<span> sortingCount (_banks, _options);

	// We create a custom count processor and give it to the sorting count algorithm
	SimkaDistanceParam distanceParams(_options);
	_processor = new SimkaCountProcessor<span> (*_stats, _nbBanks, _abundanceThreshold, _solidKind, _soliditySingle, _progress, distanceParams);
	_processor->use();
	sortingCount.addProcessor (_processor);

	// We launch the algorithm
	sortingCount.execute();


}

template<size_t span>
void SimkaAlgorithm<span>::outputMatrix(){
	_stats->outputMatrix(_outputDir, _bankNames);
}

/*
template<size_t span>
void SimkaAlgorithm<span>::outputHeatmap(){
	cout << endl << endl;
	__outputHeatmap("heatmap_presenceAbsence_sorensen", "mat_presenceAbsence_sorensen", "mat_presenceAbsence_sorensen");
	__outputHeatmap("heatmap_presenceAbsence_jaccard", "mat_presenceAbsence_jaccard_asym", "mat_presenceAbsence_jaccard");
	__outputHeatmap("heatmap_abundance_jaccard", "mat_abundance_jaccard_asym", "mat_abundance_jaccard");
	__outputHeatmap("heatmap_abundance_brayCurtis", "mat_abundance_brayCurtis", "mat_abundance_brayCurtis");
}


template<size_t span>
void SimkaAlgorithm<span>::__outputHeatmap(const string& outputFilenamePrefix, const string& matrixAsymFilename, const string& matrixNormFilename){



	string asymFilename = matrixAsymFilename + _outputFilenameSuffix + ".csv";
	string normFilename = matrixNormFilename + _outputFilenameSuffix + ".csv";
	string outputFilename = outputFilenamePrefix + _outputFilenameSuffix + ".png";


	string command = "Rscript ./Rscripts/heatmap.r " + _outputDir + "/" + asymFilename + " " + _outputDir + "/" + normFilename + " " + _outputDir + "/" + outputFilename;
	//cout << command << endl;

    try
    {
    	system(command.c_str());
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        //return EXIT_FAILURE;
    }

    //if(linePartList[0] == "dks")
    //	_heatmapDksFilename = outputFilename;
    //else
    //	_heatmapAksFilename = outputFilename;

}*/


template<size_t span>
void SimkaAlgorithm<span>::print(){

	cout << "Output folder:   " << _outputDir << endl;

	/*
	cout << "Similarity matrix:" << endl;
	cout << "\t" << "DKS (presence/absence)" << endl;
	cout << "\t\t" << "asym: " << _outputDir + "/" + _matDksPercFilename << endl;
	cout << "\t\t" << "norm: " << _outputDir + "/" + _matDksNormFilename << endl;
	cout << "\t" << "AKS (abundance)" << endl;
	cout << "\t\t" << "asym: " << _outputDir + "/" + _matAksPercFilename << endl;
	cout << "\t\t" << "norm: " << _outputDir + "/" + _matAksNormFilename << endl;

	cout << "Heatmaps:" << endl;
	cout << "\t" << "DKS (presence/absence):" << _outputDir + "/" + _heatmapDksFilename << endl;
	cout << "\t" << "AKS (abundance):" << _outputDir + "/" + _heatmapAksFilename << endl;*/


}


template<size_t span>
void SimkaAlgorithm<span>::clear(){

	if(_banks){
		//_banks->finalize();
		//delete _banks;
	}

	System::file().remove(_banksInputFilename);
    if(_processor) _processor->forget();

    for(size_t i=0; i<_tempFilenamesToDelete.size(); i++){
    	System::file().remove(_tempFilenamesToDelete[i]);
    }

    if(_stats) delete _stats;
    //if(_simkaDistance) delete _simkaDistance;
	//_banks->remove();
	//delete _processor;
}
