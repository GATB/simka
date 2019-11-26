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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_





#include "SimkaMinCommons.hpp"
#include <mutex>

class KmerSpectrumIterator{
public:

	//FILE * _is;
	//ifstream _kmerSpectrumFile;
	size_t _sketchSize;
	bool _isDone;
	u_int64_t _nbItems;

	size_t _datasetId;
	vector<vector<KmerAndCountType> >& _kmercountSketches;

	size_t _nbDatasetOffset;
	//u_int8_t* _buffer;
	//u_int64_t _bufferSize;

	//KmerAndCountType* _buffer;

	//vector<vector<KmerAndCountType> > _buffers;
	//vector<bool> _buffers_isInit;

	KmerSpectrumIterator(const string& filename, vector<vector<KmerAndCountType> >& kmercountSketches, size_t nbDatasetOffset) :
	_kmercountSketches(kmercountSketches)
	{

		//_buffer = 0;
		//_is = fopen((filename).c_str(), "rb");

		//_kmerSpectrumFile.open(filename + ".kmers", ios::binary);
		//_sketchSize = sketchSize;
		_nbDatasetOffset = nbDatasetOffset;
		//cout << nbDatasetOffset << endl;

        //_buffer  = (KmerAndCountType*) MALLOC (sizeof(KmerAndCountType) * _sketchSize);
	}

	~KmerSpectrumIterator(){
		//fclose(_is);
		//_kmerSpectrumFile.close();
		//if(_buffer){FREE (_buffer);}
	}

	void first(size_t datasetId){
		//cout << datasetId << " " << _nbDatasetOffset << endl;
		_datasetId = datasetId-_nbDatasetOffset;
		//cout << datasetId << " " << _nbDatasetOffset << " " << _datasetId<< endl;
		//if(_buffer){FREE (_buffer);}
		//u_int64_t pos = KMER_SPECTRUM_HEADER_SIZE + (datasetId*_sketchSize*sizeof(KmerAndCountType));
		//fseek(_is, pos, SEEK_SET);
		_nbItems = 0;
		_sketchSize =  _kmercountSketches[_datasetId].size();
		//cout << sizeof(KmerAndCountType) << endl;
		//_kmerSpectrumFile.read((char*)_buffer, 10*_sketchSize);
		//int res = fread(_buffer, sizeof(KmerAndCountType), _sketchSize, _is);
	}

	inline bool isDone(){
		return _nbItems >= _sketchSize;
	}

	inline void next(u_int64_t& kmer, KmerCountType& count){
		//KmerAndCountType kmerCount; // =  _buffer[_nbItems];
		//memcpy(&kmerCount, &_buffer[_nbItems*10], 10);

		KmerAndCountType& kmerCount = _kmercountSketches[_datasetId][_nbItems];
		//cout << _datasetId << " " << _nbDatasetOffset << " " << _kmercountSketches[_datasetId].size() << endl;
		kmer = kmerCount._kmer;
		count = kmerCount._count;

		//cout << kmer << " " << count << endl;

		//_kmerSpectrumFile.read((char*)(&kmer), sizeof(kmer));
		//_kmerSpectrumFile.read((char*)(&count), sizeof(count));
		_nbItems += 1;
	}

};

class ComputeDistanceManager{
public:

	ofstream& _distanceMatrixJaccard;
	ofstream& _distanceMatrixBrayCurtis;

	//size_t _sketchSize;
	KmerSpectrumIterator* _kmerSpectrumiterator1;
	KmerSpectrumIterator* _kmerSpectrumiterator2;
	bool _isSymmetrical;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbDistinctSharedKmers;
	u_int64_t _nbKmers;
	u_int64_t _nbSharedKmers;

	size_t _nbDatasets1;
	size_t _nbDatasets2;

	vector<PairwiseDistance> _jaccardDistances;
	vector<PairwiseDistance> _braycurtisDistances;
	u_int64_t _jaccardDistances_nb;
	mutex& _mutex;


	//u_int64_t nbLala;

	ComputeDistanceManager(const string& filename1, const string& filename2, ofstream& distanceMatrixJaccard, ofstream& distanceMatrixBrayCurtis, bool isSymmetrical, size_t nbDatasets1, size_t nbDatasets2, mutex& mutex, size_t main_start_i, size_t main_start_j, size_t n_i, size_t n_j, vector<vector<KmerAndCountType> >& _kmercountSketches_i, vector<vector<KmerAndCountType> >& _kmercountSketches_j)
	: _distanceMatrixJaccard(distanceMatrixJaccard), _distanceMatrixBrayCurtis(distanceMatrixBrayCurtis), _mutex(mutex)
	{
		//_sketchSize = sketchSize;
		_kmerSpectrumiterator1 = new KmerSpectrumIterator(filename1, _kmercountSketches_i, main_start_i);
		_kmerSpectrumiterator2 = new KmerSpectrumIterator(filename2, _kmercountSketches_j, main_start_j);
		_isSymmetrical = isSymmetrical;
		_nbDatasets1 = nbDatasets1;
		_nbDatasets2 = nbDatasets2;

		_jaccardDistances.resize(1000);
		_braycurtisDistances.resize(1000);
		_jaccardDistances_nb = 0;

		//nbLala = 0;
	}

	~ComputeDistanceManager(){
		delete _kmerSpectrumiterator1;
		delete _kmerSpectrumiterator2;

		if(_jaccardDistances_nb > 0){
			writeDistances();
			/*
			_mutex.lock();
			for(size_t i=0; i<_jaccardDistances.size() ; i++){
				PairwiseDistance& jaccard = _jaccardDistances[i];
				PairwiseDistance& braycurtis = _braycurtisDistances[i];

				u_int64_t pos = jaccard._i*_nbDatasets2*sizeof(DistanceValueType) + (jaccard._j*sizeof(DistanceValueType));
				_distanceMatrixJaccard.seekp(pos);
				_distanceMatrixBrayCurtis.seekp(pos);
				_distanceMatrixJaccard.write((const char*)&jaccard._distance, sizeof(jaccard._distance));
				_distanceMatrixBrayCurtis.write((const char*)&braycurtis._distance, sizeof(braycurtis._distance));

				if(_isSymmetrical){
					u_int64_t pos = jaccard._j*_nbDatasets1*sizeof(DistanceValueType) + (jaccard._i*sizeof(DistanceValueType));
					_distanceMatrixJaccard.seekp(pos);
					_distanceMatrixBrayCurtis.seekp(pos);
					_distanceMatrixJaccard.write((const char*)&jaccard._distance, sizeof(jaccard._distance));
					_distanceMatrixBrayCurtis.write((const char*)&braycurtis._distance, sizeof(braycurtis._distance));
				}
			}
			_mutex.unlock();
			_jaccardDistances_nb = 0;
			//_braycurtisDistances.clear();
			//_jaccardDistances.clear();

			 */
		}

		//cout << nbLala << endl;
	}

	void computeDistance_unsynch(size_t i, size_t j){

		//nbLala += 1;
		//_mutex.lock();
		//lala += 1;
		//_mutex.unlock();

		_nbDistinctSharedKmers = 0;
		_nbDistinctKmers = 0;
		_nbKmers = 0;
		_nbSharedKmers = 0;

		_kmerSpectrumiterator1->first(i);
		_kmerSpectrumiterator2->first(j);

		u_int64_t sketchSize = min(_kmerSpectrumiterator1->_sketchSize, _kmerSpectrumiterator2->_sketchSize);

		u_int64_t kmer1;
		u_int64_t kmer2;
		KmerCountType count1;
		KmerCountType count2;

		_kmerSpectrumiterator1->next(kmer1, count1);
		_kmerSpectrumiterator2->next(kmer2, count2);

		while(_nbDistinctKmers < sketchSize){ //_nbDistinctKmers < _sketchSize && (!_kmerSpectrumiterator1->isDone()) && (!_kmerSpectrumiterator2->isDone()) ){
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

		DistanceValueType jaccard;
		DistanceValueType braycurtis;


		if(_nbDistinctKmers == 0){
			jaccard = 1;
		}
		else{
			jaccard = 1 - (long double) _nbDistinctSharedKmers / (long double) _nbDistinctKmers;
		}

		if(_nbKmers == 0){
			braycurtis = 1;
		}else{
			braycurtis =  1 - (long double) (2*_nbSharedKmers) / (long double) _nbKmers;
		}


		//_mutex.lock();
		//cout << i << " " << j << " " << braycurtis << endl;
		//_mutex.unlock();

		_jaccardDistances[_jaccardDistances_nb].set(i, j, jaccard);
		_braycurtisDistances[_jaccardDistances_nb].set(i, j, braycurtis);
		_jaccardDistances_nb += 1;

		if(_jaccardDistances_nb == _jaccardDistances.size()){
			writeDistances();
			//_braycurtisDistances.clear();
			//_jaccardDistances.clear();
		}

		//cout << "NB DISTINCT KMERS: " << _nbDistinctKmers << endl;
		//cout << "NB SHARED DISTINCT KMERS: " << _nbDistinctSharedKmers << endl;
		//cout << "JACCARD: " <<  << endl;
		//cout << "BRAY CURTIS: " << 1 - (long double) (2*_nbSharedKmers) / (long double) _nbKmers << endl;
	}

	void writeDistances(){
		_mutex.lock();
		for(size_t i=0; i<_jaccardDistances_nb ; i++){
			PairwiseDistance& jaccard = _jaccardDistances[i];
			PairwiseDistance& braycurtis = _braycurtisDistances[i];

			u_int64_t pos = jaccard._i*_nbDatasets2*sizeof(DistanceValueType) + (jaccard._j*sizeof(DistanceValueType));
			_distanceMatrixJaccard.seekp(pos);
			_distanceMatrixBrayCurtis.seekp(pos);
			_distanceMatrixJaccard.write((const char*)&jaccard._distance, sizeof(jaccard._distance));
			_distanceMatrixBrayCurtis.write((const char*)&braycurtis._distance, sizeof(braycurtis._distance));

			if(_isSymmetrical){
				u_int64_t pos = jaccard._j*_nbDatasets1*sizeof(DistanceValueType) + (jaccard._i*sizeof(DistanceValueType));
				_distanceMatrixJaccard.seekp(pos);
				_distanceMatrixBrayCurtis.seekp(pos);
				_distanceMatrixJaccard.write((const char*)&jaccard._distance, sizeof(jaccard._distance));
				_distanceMatrixBrayCurtis.write((const char*)&braycurtis._distance, sizeof(braycurtis._distance));
			}
		}
		_mutex.unlock();
		_jaccardDistances_nb = 0;
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

	u_int32_t _sketchSize_1, _sketchSize_2;
	u_int32_t _seed;
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
	u_int32_t _nbDataset1;
	u_int32_t _nbDataset2;

	ofstream _distanceMatrixJaccard;
	ofstream _distanceMatrixBrayCurtis;

	mutex _mutex;

	IteratorListener* _progress;
	u_int64_t _progress_distanceStep;
	//u_int64_t _progress_nbDistancesToCompute;
	u_int64_t _progress_nbDistancesComputed;
	//string _progress_text;

	size_t _start_i, _start_j;
	size_t _n_i, _n_j;

	vector<vector<KmerAndCountType> > _kmercountSketches_i;
	vector<vector<KmerAndCountType> > _kmercountSketches_j;


	SimkaMinDistanceAlgorithm(IProperties* options):
		Algorithm("simkaMinDistanceAlgorithm", -1, options)
	{
	}

	void execute(){
		//pthread_mutex_init(&_mutex, NULL);

		parseArgs();

		readInfos();

		loadSketches();
		distance();
		//createDirs();
		//SimkaCommons::checkInputValidity(_outputDirTemp, _inputFilename);
		//countDatasets();


		//string command = "rm -rf " + _outputDirTemp;
		//system(command.c_str());

		cout << "Output results: " << _outputDir << endl;
	}

	void parseArgs(){

		_options = getInput();

		//_sketchSize = _options->getInt(STR_SIMKA_SKETCH_SIZE);
	    _nbCores = _options->getInt(STR_NB_CORES);
		_inputFilename1 = _options->getStr(STR_SIMKA_URI_INPUT_1);
		_inputFilename2 = _options->getStr(STR_SIMKA_URI_INPUT_2);
		_outputDir = _options->getStr(STR_URI_OUTPUT);

		_start_i = _options->getInt("-start-i");
		_start_j = _options->getInt("-start-j");
		_n_i = _options->getInt("-n-i");
		_n_j = _options->getInt("-n-j");

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
		//_nbDataset1 = SimkaMinCommons::readNbDatasets(_inputFilename1);
		//_nbDataset2 = SimkaMinCommons::readNbDatasets(_inputFilename2);

		//u_int32_t sketchSize1;
		//u_int32_t sketchSize2;
		u_int8_t kmerSizeDummy;
		SimkaMinCommons::getKmerInfos(_inputFilename1, kmerSizeDummy, _sketchSize_1, _seed, _nbDataset1);
		SimkaMinCommons::getKmerInfos(_inputFilename2, kmerSizeDummy, _sketchSize_2, _seed, _nbDataset2);
		//_sketchSize = min(sketchSize1, sketchSize2);



		if(_sketchSize_1 != _sketchSize_2){
			cout << "WARNING: both spectrums have different sizes (" << _sketchSize_1 << " and " << _sketchSize_2 << "), will use " << min(_sketchSize_1, _sketchSize_2)  << " k-mers" << endl;
		}

		if(_n_i == 0){
			_n_i = _nbDataset1;
		}
		if(_n_j == 0){
			_n_j = _nbDataset2;
		}
		//_nbdatasetsToProcess = min(_nbdatasetsToProcess, )
		//cout << _nbDataset1 << " " << _nbDataset2 << endl;
		//cout << _sketchSize << endl;
		//cout << _seed << endl;
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

	void loadSketches(){

		ifstream sketchFile_1;
		sketchFile_1.open(_inputFilename1.c_str(), ios::binary);
		ifstream sketchFile_2;
		sketchFile_2.open(_inputFilename2.c_str(), ios::binary);

		_kmercountSketches_i.resize(_n_i);
		_kmercountSketches_j.resize(_n_j);

		u_int32_t sketchSize = min(_sketchSize_1, _sketchSize_2);

		size_t index = 0;
		for(size_t i=_start_i; i<_start_i+_n_i; i++){

			u_int64_t pos = KMER_SPECTRUM_HEADER_SIZE + (i*_sketchSize_1*sizeof(KmerAndCountType));
			sketchFile_1.seekg(pos);

			_kmercountSketches_i[index].resize(sketchSize);
			//for(size_t k=0; k<_sketchSize; k++){
				sketchFile_1.read((char*)&(_kmercountSketches_i[index][0]), sizeof(KmerAndCountType)*sketchSize);
			//}
			index += 1;


		}

		index = 0;
		for(size_t j=_start_j; j<_start_j+_n_j; j++){

			u_int64_t pos = KMER_SPECTRUM_HEADER_SIZE + (j*_sketchSize_2*sizeof(KmerAndCountType));
			sketchFile_2.seekg(pos);

			_kmercountSketches_j[index].resize(sketchSize);
			sketchFile_2.read((char*)&(_kmercountSketches_j[index][0]), sizeof(KmerAndCountType)*sketchSize);
			//for(size_t k=0; k<_sketchSize; k++){
				//sketchFile_2.read(&_kmercountSketches_j[index][k], sizeof(KmerAndCountType));
			//}
			index += 1;

		}

		sketchFile_1.close();
		sketchFile_2.close();


		for(size_t i=0; i<_kmercountSketches_i.size(); i++){
			u_int64_t start = 0;
			for(size_t j=0; j<_sketchSize_1; j++){
				if(_kmercountSketches_i[i][j]._kmer == 0){
					start += 1;
				}
			}
			_kmercountSketches_i[i].erase(_kmercountSketches_i[i].begin(), _kmercountSketches_i[i].begin()+start);
		}

		for(size_t i=0; i<_kmercountSketches_j.size(); i++){
			u_int64_t start = 0;
			for(size_t j=0; j<_kmercountSketches_j[i].size(); j++){
				if(_kmercountSketches_j[i][j]._kmer == 0){
					start += 1;
				}
			}
			_kmercountSketches_j[i].erase(_kmercountSketches_j[i].begin(), _kmercountSketches_j[i].begin()+start);
		}
	}

	void distance(){

		if(System::file().doesExist(_outputDir + "/mat_presenceAbsence_jaccard.bin")){
			_distanceMatrixJaccard.open((_outputDir + "/mat_presenceAbsence_jaccard.bin").c_str(), ios::binary | ios::in);
			_distanceMatrixBrayCurtis.open((_outputDir + "/mat_abundance_braycurtis.bin").c_str(), ios::binary | ios::in);
		}
		else{
			_distanceMatrixJaccard.open((_outputDir + "/mat_presenceAbsence_jaccard.bin").c_str(), ios::binary);
			_distanceMatrixBrayCurtis.open((_outputDir + "/mat_abundance_braycurtis.bin").c_str(), ios::binary);
		}

		bool isSymmetrical = false;
		if(_inputFilename1 == _inputFilename2 && _start_i == _start_j){
			computeDistanceSymetrical();
			isSymmetrical = true;
		}
		else{
			computeDistanceRectangle();
		}


		for(size_t i=0; i<_threads.size(); i++){
			_threads[i]->join();
			delete _threads[i];
			//cout << i << endl;
		}

		_progress->finish();

		//Fill diagonal with 0

		if(isSymmetrical){
			for(size_t i=_start_i; i<_start_i+_n_i; i++){
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

		//string command = "cp " + string(_inputFilename1+".ids") + " " + _outputDir + "/matrix_infos.ids ";
		//cout << command << endl;
		//system(command.c_str());
	}



	void computeDistanceSymetrical(){
		//cout << "compute symetrical distances" << endl;

		u_int64_t nbDistancesToCompute = (_n_i*(_n_i-1)) / 2;
		//u_int64_t nbDistancesToCompute = _nbDataset1*_nbDataset1; //(_nbDataset1*(_nbDataset1-1)) / 2;
		u_int64_t nbDistancePerThreads = nbDistancesToCompute / _nbCores;
		u_int64_t nbDistancesRemaining = nbDistancesToCompute-(nbDistancePerThreads*_nbCores);

		//vector<size_t> startDistanceI;
		//vector<size_t> startDistanceJ;
		//size_t si=0;
		//size_t sj=0;

		//cout << "NB CORES: " << _nbCores << endl;
		//cout << "NB DISTANCES: " << nbDistancesToCompute << endl;
		//cout << "NB DISTANCES PER CORE: " << nbDistancePerThreads << endl;

		_progress = this->createIteratorListener (nbDistancesToCompute, "Computing distances");
		_progress->init ();
		_progress_distanceStep = max((u_int64_t)1, (u_int64_t) (nbDistancePerThreads / 100));

		u_int64_t nbDistances = 0;

		size_t nbRunnedThreads = 0;
		size_t i=_start_i;
		size_t j=i+1;
		size_t maxDatasets = _start_i+_n_i;//min((u_int64_t)_start_i+_nbdatasetsToProcess, (u_int64_t)_nbDataset1);


		//_computeDistanceManagers.push_back();
		thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_unsynch, this, i, j, nbDistancePerThreads, nbRunnedThreads);
		_threads.push_back(t);
		//computeDistances_unsynch(i, j, nbDistancePerThreads, true);

		bool done = false;
		nbRunnedThreads += 1;

		for(; i<maxDatasets; i++){

			if(done) break;

			for(j=i+1; j<maxDatasets; j++){

				if(done) break;

				if(nbDistances >= nbDistancePerThreads){
					//cout << i << " " << j << endl;

					//cout << "lol: " << nbRunnedThreads << " " << nbDistancesRemaining << endl;
					if(nbRunnedThreads == _nbCores-1){ //Last threads compute remaining distances
						//cout << " LOL " << endl;);
						thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_unsynch, this, i, j, nbDistancePerThreads+nbDistancesRemaining, nbRunnedThreads);
						_threads.push_back(t);
						//computeDistances_unsynch(i, j, nbDistancePerThreads+nbDistancesRemaining, true);
						done = true;
						//nbDistances -= nbDistancesRemaining;
					}
					else{
						thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_unsynch, this, i, j, nbDistancePerThreads, nbRunnedThreads);
						_threads.push_back(t);
						//computeDistances_unsynch(i, j, nbDistancePerThreads, true);
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

		//cout << "compute rectangle distances" << endl;


		u_int64_t nbDistancesToCompute = _n_i*_n_j;
		//u_int64_t nbDistancesToCompute = _nbDataset1*_nbDataset1; //(_nbDataset1*(_nbDataset1-1)) / 2;
		u_int64_t nbDistancePerThreads = nbDistancesToCompute / _nbCores;
		u_int64_t nbDistancesRemaining = nbDistancesToCompute-(nbDistancePerThreads*_nbCores);

		//vector<size_t> startDistanceI;
		//vector<size_t> startDistanceJ;
		//size_t si=0;
		//size_t sj=0;

		//cout << "NB CORES: " << _nbCores << endl;
		//cout << "NB DISTANCES: " << nbDistancesToCompute << endl;
		//cout << "NB DISTANCES PER CORE: " << nbDistancePerThreads << endl;
		//cout << "NB DISTANCES REMAINING: " << nbDistancesRemaining << endl;

		_progress = this->createIteratorListener (nbDistancesToCompute, "Computing distances");
		_progress->init ();

		u_int64_t nbDistances = 0;

		size_t nbRunnedThreads = 0;
		size_t i=_start_i;
		size_t j=_start_j;


		thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_rectanglular_unsynch, this, i, j, nbDistancePerThreads, nbRunnedThreads);
		_threads.push_back(t);

		bool done = false;
		nbRunnedThreads += 1;

		//size_t maxDatasetsI = min((u_int64_t)_start_i+_nbdatasetsToProcess, (u_int64_t)_nbDataset1);
		for(i=_start_i; i<_start_i+_n_i; i++){

			if(done) break;

			//size_t maxDatasetsJ = min((u_int64_t)_start_j+_nbdatasetsToProcess, (u_int64_t)_nbDataset2);
			for(j=_start_j; j<_start_j+_n_j; j++){

				if(done) break;

				if(nbDistances >= nbDistancePerThreads){
					//cout << i << " " << j << endl;

					//cout << "lol: " << nbRunnedThreads << " " << nbDistancesRemaining << endl;
					if(nbRunnedThreads == _nbCores-1){ //Last threads compute remaining distances
						//cout << " LOL " << endl;);
						//cout << i << " " << j << endl;
						thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_rectanglular_unsynch, this, i, j, nbDistancePerThreads+nbDistancesRemaining, nbRunnedThreads);
						_threads.push_back(t);
						//computeDistances_unsynch(i, j, nbDistancePerThreads+nbDistancesRemaining, true);
						done = true;
						//nbDistances -= nbDistancesRemaining;
					}
					else{
						//cout << i << " " << j << endl;
						thread* t = new thread(&SimkaMinDistanceAlgorithm::computeDistances_rectanglular_unsynch, this, i, j, nbDistancePerThreads, nbRunnedThreads);
						_threads.push_back(t);
						//computeDistances_unsynch(i, j, nbDistancePerThreads, true);
					}

					nbRunnedThreads += 1;
					nbDistances = 0;
				}
				//cout << nbDistances << " " << nbDistancePerThreads << endl;
				nbDistances += 1;
			}
		}

	}

	void computeDistances_unsynch(size_t si, size_t sj, size_t nbDistancesToCompute, size_t id){

		ComputeDistanceManager computeDistanceManager(_inputFilename1, _inputFilename2, _distanceMatrixJaccard, _distanceMatrixBrayCurtis, true, _nbDataset1, _nbDataset2, _mutex, _start_i, _start_j, _n_i, _n_j, _kmercountSketches_i, _kmercountSketches_j);

		//cout << "-------------------" << endl;
		u_int64_t progress_nbComputedistances = 0;
		u_int64_t nbComputedDistances = 0;


		//size_t maxDatasetsI = min((u_int64_t)si+_nbdatasetsToProcess, (u_int64_t)_nbDataset1);
		for(size_t j=sj; j<_start_i+_n_i; j++){
			//cout << j << endl;
			//cout << si << " " << j << endl;
			computeDistanceManager.computeDistance_unsynch(si, j);
			nbComputedDistances += 1;
			progress_nbComputedistances += 1;
			if(nbComputedDistances >= nbDistancesToCompute) break;
		}
		si += 1;

		if(nbComputedDistances < nbDistancesToCompute){

			//cout << "lala2" << endl;
			for(size_t i=si; i<_start_i+_n_i; i++){
				for(size_t j=i+1; j<_start_i+_n_i; j++){
					//cout << i << " " << j << endl;
					computeDistanceManager.computeDistance_unsynch(i, j);
					nbComputedDistances += 1;
					progress_nbComputedistances += 1;

					//_mutex.lock();
					//cout << progress_nbComputedistances << "   " << _progress_distanceStep << endl;
					//_mutex.unlock();

					if(progress_nbComputedistances > _progress_distanceStep){
						_mutex.lock();
						_progress->inc(progress_nbComputedistances);
						_mutex.unlock();
						progress_nbComputedistances = 0;
					}
					if(nbComputedDistances >= nbDistancesToCompute) break;
				}
				if(nbComputedDistances >= nbDistancesToCompute) break;
			}
		}


		_mutex.lock();
		_progress->inc(progress_nbComputedistances);
		_mutex.unlock();

	}

	void computeDistances_rectanglular_unsynch(size_t si, size_t sj, size_t nbDistancesToCompute, size_t id){

		//isSymetrical set to false
		ComputeDistanceManager computeDistanceManager(_inputFilename1, _inputFilename2, _distanceMatrixJaccard, _distanceMatrixBrayCurtis, false, _nbDataset1, _nbDataset2, _mutex, _start_i, _start_j, _n_i, _n_j, _kmercountSketches_i, _kmercountSketches_j);

		//_mutex.lock();
		//cout << "------------------- " << si << " " << sj << endl;
		//_mutex.unlock();

		u_int64_t nbComputedDistances = 0;
		u_int64_t progress_nbComputedistances = 0;


		//size_t maxDatasetsI = min((u_int64_t)si+_nbdatasetsToProcess, (u_int64_t)_nbDataset1);
		//size_t maxDatasetsJ = min((u_int64_t)sj+_nbdatasetsToProcess, (u_int64_t)_nbDataset2);
		for(size_t j=sj; j<_start_j+_n_j; j++){
			//cout << si << " " << j << endl;
			computeDistanceManager.computeDistance_unsynch(si, j);
			nbComputedDistances += 1;
			progress_nbComputedistances += 1;
			if(nbComputedDistances >= nbDistancesToCompute) break;
		}
		si += 1;

		if(nbComputedDistances < nbDistancesToCompute){

			for(size_t i=si; i<_start_i+_n_i; i++){
				for(size_t j=_start_j; j<_start_j+_n_j; j++){ // (0 instead of i+1)
					//cout << i << " " << j << endl;
					computeDistanceManager.computeDistance_unsynch(i, j);
					nbComputedDistances += 1;
					progress_nbComputedistances += 1;

					if(progress_nbComputedistances > _progress_distanceStep){
						_mutex.lock();
						_progress->inc(progress_nbComputedistances);
						_mutex.unlock();
						progress_nbComputedistances = 0;
					}
					if(nbComputedDistances >= nbDistancesToCompute) break;
				}
				if(nbComputedDistances >= nbDistancesToCompute) break;
			}
		}


		_mutex.lock();
		//cout << nbComputedDistances << " " << nbDistancesToCompute << endl;
		_progress->inc(progress_nbComputedistances);
		_mutex.unlock();

	}

};

















class SimkaMinDistance : public Tool{
public:


	SimkaMinDistance(): Tool ("SimkaMin-Distance"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for distance matrices", false, "./simkaMin_results"));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_2, "filename to a sketch file to compare with -in1", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_1, "filename to a sketch file to compare with -in2", true));

	    parser->push_back (new OptionOneParam ("-start-i", "start i (row)", false, "0"));
	    parser->push_back (new OptionOneParam ("-start-j", "start j (column)", false, "0"));
	    parser->push_back (new OptionOneParam ("-n-i", "Nb datasets to process (row)", false, "0"));
	    parser->push_back (new OptionOneParam ("-n-j", "Nb datasets to process (column)", false, "0"));

	}


	void execute ()
	{
		IProperties* args = getInput();

		u_int32_t seed1;
		u_int32_t seed2;
		u_int32_t dummy;
		u_int8_t kmerSize1;
		u_int8_t kmerSize2;
		string inputFilename1 = args->getStr(STR_SIMKA_URI_INPUT_1);
		string inputFilename2 = args->getStr(STR_SIMKA_URI_INPUT_2);
		SimkaMinCommons::getKmerInfos(inputFilename1, kmerSize1, dummy, seed1, dummy);
		SimkaMinCommons::getKmerInfos(inputFilename2, kmerSize2, dummy, seed2, dummy);
		//size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

		if(kmerSize1 != kmerSize2){
			cerr << "ERROR: can't compare both sketches because of different kmer sizes (" << kmerSize1 << " and " << kmerSize2 << ")" << endl;
			exit(1);
		}

		if(seed1 != seed2){
			cerr << "ERROR: can't compare both sketches because of different seeds (" << seed1 << " and " << seed2 << ")" << endl;
			exit(1);
		}

		//cout << seed1 << " " << seed2 << endl;
		SimkaMinDistanceAlgorithm* algo = new SimkaMinDistanceAlgorithm(args);
		algo->execute();
		delete algo;
	}

};




#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_ */
