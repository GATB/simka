/*
 * SimkaMinDistance.hpp
 *
 *  Created on: 27 juin 2017
 *      Author: gbenoit
 */

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_





#include "SimkaMinCommons.hpp"


class KmerSpectrumIterator{
public:

	FILE * _is;
	//ifstream _kmerSpectrumFile;
	size_t _sketchSize;
	bool _isDone;
	u_int64_t _nbItems;

	//u_int8_t* _buffer;
	//u_int64_t _bufferSize;

	KmerAndCountType* _buffer;

	KmerSpectrumIterator(const string& filename, size_t sketchSize){
		_buffer = 0;
		_is = fopen((filename + ".kmers").c_str(), "rb");

		//_kmerSpectrumFile.open(filename + ".kmers", ios::binary);
		_sketchSize = sketchSize;


        _buffer  = (KmerAndCountType*) MALLOC (sizeof(KmerAndCountType) * _sketchSize);
		//_bufferSize = _sketchSize * (sizeof(KmerAndCountType));
		//_buffer = new u_int8_t[_bufferSize];
		//_buffer.resize(_bufferSize);
		//cout << _sketchSize << " " << sizeof(KmerAndCountType) << " " << _bufferSize << endl;
	}

	~KmerSpectrumIterator(){
		fclose(_is);
		//_kmerSpectrumFile.close();
		if(_buffer){FREE (_buffer);}
	}

	void first(size_t datasetId){
		//if(_buffer){FREE (_buffer);}
		u_int64_t pos = KMER_SPECTRUM_HEADER_SIZE + (datasetId*_sketchSize*(sizeof(u_int64_t)+sizeof(KmerCountType)));
		fseek(_is, pos, SEEK_SET);
		_nbItems = 0;
		//cout << sizeof(KmerAndCountType) << endl;
		//_kmerSpectrumFile.read((char*)_buffer, 10*_sketchSize);
		int res = fread(_buffer, sizeof(KmerAndCountType), _sketchSize, _is);
	}

	inline bool isDone(){
		return _nbItems >= _sketchSize;
	}

	inline void next(u_int64_t& kmer, KmerCountType& count){
		//KmerAndCountType kmerCount; // =  _buffer[_nbItems];
		//memcpy(&kmerCount, &_buffer[_nbItems*10], 10);

		KmerAndCountType kmerCount = _buffer[_nbItems];
		kmer = kmerCount._kmer;
		count = kmerCount._count;
		cout << kmer << " " << count << endl;

		//_kmerSpectrumFile.read((char*)(&kmer), sizeof(kmer));
		//_kmerSpectrumFile.read((char*)(&count), sizeof(count));
		_nbItems += 1;
	}

};

class ComputeDistanceManager{
public:

	ofstream& _distanceMatrixJaccard;
	ofstream& _distanceMatrixBrayCurtis;

	size_t _sketchSize;
	KmerSpectrumIterator* _kmerSpectrumiterator1;
	KmerSpectrumIterator* _kmerSpectrumiterator2;
	bool _isSymmetrical;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbDistinctSharedKmers;
	u_int64_t _nbKmers;
	u_int64_t _nbSharedKmers;

	size_t _nbDatasets1;
	size_t _nbDatasets2;

	ComputeDistanceManager(const string& filename1, const string& filename2, size_t sketchSize, ofstream& distanceMatrixJaccard, ofstream& distanceMatrixBrayCurtis, bool isSymmetrical, size_t nbDatasets1, size_t nbDatasets2)
	: _distanceMatrixJaccard(distanceMatrixJaccard), _distanceMatrixBrayCurtis(distanceMatrixBrayCurtis)
	{
		_sketchSize = sketchSize;
		_kmerSpectrumiterator1 = new KmerSpectrumIterator(filename1, _sketchSize);
		_kmerSpectrumiterator2 = new KmerSpectrumIterator(filename2, _sketchSize);
		_isSymmetrical = isSymmetrical;
		_nbDatasets1 = nbDatasets1;
		_nbDatasets2 = nbDatasets2;
	}

	~ComputeDistanceManager(){
		delete _kmerSpectrumiterator1;
		delete _kmerSpectrumiterator2;
	}

	void computeDistance_unsynch(size_t i, size_t j){

		_nbDistinctSharedKmers = 0;
		_nbDistinctKmers = 0;
		_nbKmers = 0;
		_nbSharedKmers = 0;

		_kmerSpectrumiterator1->first(i);
		_kmerSpectrumiterator2->first(j);

		u_int64_t kmer1;
		u_int64_t kmer2;
		KmerCountType count1;
		KmerCountType count2;

		_kmerSpectrumiterator1->next(kmer1, count1);
		_kmerSpectrumiterator2->next(kmer2, count2);

		while(_nbDistinctKmers < _sketchSize){ //_nbDistinctKmers < _sketchSize && (!_kmerSpectrumiterator1->isDone()) && (!_kmerSpectrumiterator2->isDone()) ){
			//cout << kmer1 << " " << kmer2 << endl;
			if(kmer1 > kmer2){
				_nbDistinctKmers += 1;
				_nbKmers += count2;

				if(_kmerSpectrumiterator2->isDone()) break;
				_kmerSpectrumiterator2->next(kmer2, count2);
			}
			else if(kmer1 < kmer2){
				_nbDistinctKmers += 1;
				_nbKmers += count1;

				if(_kmerSpectrumiterator1->isDone()) break;
				_kmerSpectrumiterator1->next(kmer1, count1);
			}
			else{
				_nbDistinctKmers += 1;
				_nbKmers += count1 + count2;
				_nbDistinctSharedKmers += 1;
				_nbSharedKmers += min(count1, count2);

				if(_kmerSpectrumiterator2->isDone() || _kmerSpectrumiterator1->isDone()) break;
				_kmerSpectrumiterator1->next(kmer1, count1);
				_kmerSpectrumiterator2->next(kmer2, count2);
			}
		}

		DistanceValueType jaccard =  1 - (long double) _nbDistinctSharedKmers / (long double) _nbDistinctKmers;
		DistanceValueType braycurtis =  1 - (long double) (2*_nbSharedKmers) / (long double) _nbKmers;

		u_int64_t pos = i*_nbDatasets1*sizeof(DistanceValueType) + (j*sizeof(DistanceValueType));
		_distanceMatrixJaccard.seekp(pos);
		_distanceMatrixBrayCurtis.seekp(pos);
		_distanceMatrixJaccard.write((const char*)&jaccard, sizeof(jaccard));
		_distanceMatrixBrayCurtis.write((const char*)&braycurtis, sizeof(braycurtis));

		if(_isSymmetrical){
			u_int64_t pos = j*_nbDatasets1*sizeof(DistanceValueType) + (i*sizeof(DistanceValueType));
			_distanceMatrixJaccard.seekp(pos);
			_distanceMatrixBrayCurtis.seekp(pos);
			_distanceMatrixJaccard.write((const char*)&jaccard, sizeof(jaccard));
			_distanceMatrixBrayCurtis.write((const char*)&braycurtis, sizeof(braycurtis));
		}
		//cout << "NB DISTINCT KMERS: " << _nbDistinctKmers << endl;
		//cout << "NB SHARED DISTINCT KMERS: " << _nbDistinctSharedKmers << endl;
		//cout << "JACCARD: " <<  << endl;
		//cout << "BRAY CURTIS: " << 1 - (long double) (2*_nbSharedKmers) / (long double) _nbKmers << endl;
	}

};




class SimkaMinDistanceAlgorithm : public Algorithm
{

public:

	size_t _nbCores;
	string _outputDir;
	string _inputFilename1;
	string _inputFilename2;

	//pair<CountNumber, CountNumber> _abundanceThreshold;
	//SIMKA_SOLID_KIND _solidKind;
	//bool _soliditySingle;
	//int64_t _maxNbReads;
	//size_t _minReadSize;
	//double _minReadShannonIndex;
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


	//vector<vector<u_int64_t>> _bufferKmers;
	//vector<vector<u_int32_t>> _bufferCounts;
	//vector<size_t> _bufferIndex;

	//vector<u_int64_t> _minHashValues;
	//vector<u_int64_t> _minHashKmers;
	//vector<u_int32_t> _minHashKmersCounts;

	size_t _sketchSize;
	//pthread_mutex_t _mutex;

	//typedef typename SelectKmersCommand<span>::KmerCountSorter KmerCountSorter;
	//std::priority_queue< u_int64_t, vector<u_int64_t>, KmerCountSorter> _kmerCountSorter;
	//KmerCountDictionaryType _kmerCounts;

	size_t _nbBanks;
	//vector<string> _bankNames;
	//vector<size_t> _nbBankPerDataset;

	vector<thread*> _threads;
	size_t _maxRunningThreads;
	vector<size_t> _runningThreadIds;
	size_t _nbRunningThreads;
	vector<size_t> _finishedThreads;
	mutex countKmersMutex;

	//vector<string> _datasetIds1;
	//vector<string> _datasetIds2;
	size_t _nbDataset1;
	size_t _nbDataset2;

	ofstream _distanceMatrixJaccard;
	ofstream _distanceMatrixBrayCurtis;

	SimkaMinDistanceAlgorithm(IProperties* options):
		Algorithm("simkaMinDistanceAlgorithm", -1, options)
	{
	}

	void execute(){
		//pthread_mutex_init(&_mutex, NULL);

		parseArgs();

		readInfos();

		distance();
		//createDirs();
		//SimkaCommons::checkInputValidity(_outputDirTemp, _inputFilename);
		//countDatasets();


		//string command = "rm -rf " + _outputDirTemp;
		//system(command.c_str());

	}

	void parseArgs(){

		_options = getInput();

		//_sketchSize = _options->getInt(STR_SIMKA_SKETCH_SIZE);
	    _nbCores = _options->getInt(STR_NB_CORES);
		_inputFilename1 = _options->getStr(STR_SIMKA_URI_INPUT_1);
		_inputFilename2 = _options->getStr(STR_SIMKA_URI_INPUT_2);
		_outputDir = _options->getStr(STR_URI_OUTPUT);
		//_kmerSize = _options->getInt(STR_KMER_SIZE);


		if(!System::file().doesExist(_outputDir)){
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
				  std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
				  exit(1);
			}
		}
	}

	void readInfos(){
		_nbDataset1 = SimkaMinCommons::readNbDatasets(_inputFilename1);
		_nbDataset2 = SimkaMinCommons::readNbDatasets(_inputFilename2);

		size_t sketchSize1;
		size_t sketchSize2;
		size_t kmerSizeDummy;
		SimkaMinCommons::getKmerInfos(_inputFilename1, kmerSizeDummy, sketchSize1);
		SimkaMinCommons::getKmerInfos(_inputFilename2, kmerSizeDummy, sketchSize2);
		_sketchSize = min(sketchSize1, sketchSize2);


		if(sketchSize1 != sketchSize2){
			cout << "WARNING: both spectrums have different sizes (" << sketchSize1 << " and " << sketchSize2 << "), will use " << _sketchSize  << " k-mers" << endl;
		}
		cout << _nbDataset1 << " " << _nbDataset2 << endl;
		cout << _sketchSize << endl;
	}

	/*
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

	}*/


	void distance(){


		_distanceMatrixJaccard.open((_outputDir + "/mat_presenceAbsence_jaccard.bin").c_str(), ios::binary);
		_distanceMatrixBrayCurtis.open((_outputDir + "/mat_abundance_braycurtis.bin").c_str(), ios::binary);

		bool isSymmetrical = false;
		if(_inputFilename1 == _inputFilename2){
			computeDistanceSymetrical();
			isSymmetrical = true;
		}
		else{
			computeDistanceRectangle();
		}

		for(size_t i=0; i<_threads.size(); i++){
			_threads[i]->join();
			delete _threads[i];
		}

		//Fill diagonal with 0
		if(isSymmetrical){
			for(size_t i=0; i<_nbDataset1; i++){
				size_t j=i;
				u_int64_t pos = i*_nbDataset1*sizeof(DistanceValueType) + (j*sizeof(DistanceValueType));
				_distanceMatrixJaccard.seekp(pos);
				_distanceMatrixBrayCurtis.seekp(pos);
				DistanceValueType nullDist = 0;
				_distanceMatrixJaccard.write((const char*)&nullDist, sizeof(nullDist));
				_distanceMatrixBrayCurtis.write((const char*)&nullDist, sizeof(nullDist));
			}
		}


		_distanceMatrixJaccard.close();
		_distanceMatrixBrayCurtis.close();

		string command = "cp " + string(_inputFilename1+".ids") + " " + _outputDir + "/matrix_infos.ids ";
		cout << command << endl;
		system(command.c_str());
	}

	void computeDistanceSymetrical(){
		cout << "compute symetrical distances" << endl;

		u_int64_t nbDistancesToCompute = (_nbDataset1*(_nbDataset1-1)) / 2;
		//u_int64_t nbDistancesToCompute = _nbDataset1*_nbDataset1; //(_nbDataset1*(_nbDataset1-1)) / 2;
		u_int64_t nbDistancePerThreads = nbDistancesToCompute / _nbCores;
		u_int64_t nbDistancesRemaining = nbDistancesToCompute-(nbDistancePerThreads*_nbCores);

		//vector<size_t> startDistanceI;
		//vector<size_t> startDistanceJ;
		//size_t si=0;
		//size_t sj=0;

		cout << "NB CORES: " << _nbCores << endl;
		cout << "NB DISTANCES: " << nbDistancesToCompute << endl;
		cout << "NB DISTANCES PER CORE: " << nbDistancePerThreads << endl;


		u_int64_t nbDistances = 0;

		size_t nbRunnedThreads = 0;
		size_t i=0;
		size_t j=i+1;
		computeDistances_unsynch(i, j, nbDistancePerThreads, true);
		nbRunnedThreads += 1;
		for(; i<_nbDataset1; i++){
			for(j=i+1; j<_nbDataset1; j++){
				if(nbDistances >= nbDistancePerThreads){
					//cout << i << " " << j << endl;

					cout << "lol: " << nbRunnedThreads << " " << nbDistancesRemaining << endl;
					if(nbRunnedThreads == _nbCores){ //Last threads compute remaining distances
						cout << " LOL " << endl;
						computeDistances_unsynch(i, j, nbDistancePerThreads+nbDistancesRemaining, true);
					}
					else{
						computeDistances_unsynch(i, j, nbDistancePerThreads, true);
					}

					nbRunnedThreads += 1;
					nbDistances = 0;
				}
				nbDistances += 1;
			}
		}

		/*
		uint64_t iFloor = nbDistancePerThreads / _nbDataset1;
		uint64_t iMod = nbDistancePerThreads % _nbDataset1;

		for ( uint64_t i = 0, j = 0; i < _nbDataset1; i += iFloor, j += iMod )
		{
			if ( j >= _nbDataset1 )
			{
				if ( i == _nbDataset1 - 1 )
				{
					break;
				}

				i++;
				j -= _nbDataset1;
			}

			cout << i << " " << j << endl;


			//thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_unsynch, this, i, j);
			//_threads.push_back(t);
		}
		*/
		/*
		while(true){
			startDistanceI.push_back(si);
			startDistanceI.push_back(sj);

			u_int64_t nbDistances = 0;

			while(nbDistances < nbDistancesToCompute){
				//for(size_t i=0; i<)
			}
		}

		for(size_t i=1; i<_nbDataset1; i++){
			for(size_t j=(i+1); j<_nbDataset1; j++){

			}
		}*/
	}

	void computeDistanceRectangle(){

		cout << "compute rectangle distances" << endl;
	}

	void computeDistances_unsynch(size_t si, size_t sj, size_t nbDistancesToCompute, bool isSymmetrical){

		ComputeDistanceManager computeDistanceManager(_inputFilename1, _inputFilename2, _sketchSize, _distanceMatrixJaccard, _distanceMatrixBrayCurtis, isSymmetrical, _nbDataset1, _nbDataset2);

		cout << "-------------------" << endl;
		u_int64_t nbComputedDistances = 0;

		for(size_t j=sj; j<_nbDataset1; j++){
			cout << si << " " << j << endl;
			computeDistanceManager.computeDistance_unsynch(si, j);
			nbComputedDistances += 1;
			if(nbComputedDistances >= nbDistancesToCompute) break;
		}
		si += 1;

		if(nbComputedDistances < nbDistancesToCompute){

			for(size_t i=si; i<_nbDataset1; i++){
				for(size_t j=i+1; j<_nbDataset1; j++){
					cout << i << " " << j << endl;
					computeDistanceManager.computeDistance_unsynch(i, j);
					nbComputedDistances += 1;
					if(nbComputedDistances >= nbDistancesToCompute) break;
				}
				if(nbComputedDistances >= nbDistancesToCompute) break;
			}
		}
	}



};

















class SimkaMinDistance : public Tool{
public:


	SimkaMinDistance(): Tool ("SimkaMin-Distance"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for distance matrices", false, "./simkaMin_results"));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_2, "filename to a kmer count file to compare with -in1", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_1, "filename to a kmer count file to compare with -in2", true));

	}


	void execute ()
	{
		IProperties* args = getInput();

		size_t dummy;
		size_t kmerSize1;
		size_t kmerSize2;
		string inputFilename1 = args->getStr(STR_SIMKA_URI_INPUT_1);
		string inputFilename2 = args->getStr(STR_SIMKA_URI_INPUT_2);
		SimkaMinCommons::getKmerInfos(inputFilename1, kmerSize1, dummy);
		SimkaMinCommons::getKmerInfos(inputFilename2, kmerSize2, dummy);
		//size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

		if(kmerSize1 != kmerSize2){
			cout << "ERROR: can't compare both spectrums because of different kmer sizes (" << kmerSize1 << " and " << kmerSize2 << ")" << endl;
			exit(1);
		}

		SimkaMinDistanceAlgorithm* algo = new SimkaMinDistanceAlgorithm(args);
		algo->execute();
		delete algo;
	}

};




#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_ */
