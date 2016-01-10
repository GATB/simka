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
SimkaCountProcessor<span>::SimkaCountProcessor (SimkaStatistics& stats, size_t nbBanks, size_t kmerSize, const pair<size_t, size_t>& abundanceThreshold, SIMKA_SOLID_KIND solidKind, bool soliditySingle, double minKmerShannonIndex) :
_stats(stats)
{

	// We configure the vector for the N.(N+1)/2 possible pairs
	//_countTotal.resize (_nbBanks*(_nbBanks+1)/2);

	_nbBanks = nbBanks;
	_kmerSize = kmerSize;
	_abundanceThreshold = abundanceThreshold;
	_solidKind = solidKind;
	_soliditySingle = soliditySingle;
	_minKmerShannonIndex = minKmerShannonIndex;

	_localStats = new SimkaStatistics(_nbBanks, _stats._distanceParams);

	_nbKmerCounted = 0;
	isAbundanceThreshold = _abundanceThreshold.first > 1 || _abundanceThreshold.second < 1000000;

	_solidCounts.resize(_nbBanks);
}

template<size_t span>
SimkaCountProcessor<span>::~SimkaCountProcessor () {

	/*if(_progress){ //Simka_min
		if(_nbKmerCounted > 0){
			_progress->inc(_nbKmerCounted);
			_nbKmerCounted = 0;
		}
	}*/

	delete _localStats;
}


template<size_t span>
void SimkaCountProcessor<span>::finishClones (std::vector<ICountProcessor<span>*>& clones)
{
	cout << "finish clones" << endl;
	for (size_t i=0; i<clones.size(); i++){

		cout << i << endl;
		if (SimkaCountProcessor* clone = dynamic_cast<SimkaCountProcessor*> (clones[i])){
			cout << "cast i" << endl;
			finishClone(clone);
			//for (size_t i=0; i<this->_countTotal.size(); i++)  { this->_countTotal[i] += clone->_countTotal[i];  }
		}
	}
}

template<size_t span>
void SimkaCountProcessor<span>::finishClone(SimkaCountProcessor<span>* clone){
	cout << "finish clone" << endl;
	//cout << _stats << "   " << &clone->_stats << endl;
	_stats += *clone->_localStats;
}

template<size_t span>
bool SimkaCountProcessor<span>::isSolidVector(const CountVector& counts){

	//size_t nbBanks = 0;
	//size_t nbSolids = 0;

	for(size_t i=0; i<counts.size(); i++){

		if(counts[i] >= _abundanceThreshold.first && counts[i] <= _abundanceThreshold.second)
			return true;

		//cout << "a " <<  counts[i] << _abundanceThreshold.first << "  " << _abundanceThreshold.second << endl;
		//if(counts[i] > 0) nbBanks += 1;

	}

	//if(nbBanks > 1) return true;

	return false;

}


//template<size_t span>
//bool SimkaCountProcessor<span>::isSolid(CountNumber count){
//	return count >= _abundanceThreshold.first && count <= _abundanceThreshold.second;
//}


template<size_t span>
bool SimkaCountProcessor<span>::process (size_t partId, const Type& kmer, const CountVector& counts, CountNumber sum){




	/*
	if(_progress){ //Simka_min
		if(_nbKmerCounted > 500000){
			_progress->inc(_nbKmerCounted);
			_nbKmerCounted = 0;
		}
	}*/

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

	if(_minKmerShannonIndex != 0){
		double shannonIndex = getShannonIndex(kmer);
		if(shannonIndex < _minKmerShannonIndex){
			return false;
		}
	}

	/*
	float Ri = 500000;
	float Rtotal = Ri * _nbBanks;
	float Ntotal = _totalAbundance;
	float X2j = 0;
	for(size_t i=0; i<counts.size(); i++){

		float Ni = counts[i];

		X2j += pow((Ni/Ntotal - Ri/Rtotal), 2) / (Ri / (Rtotal*Ntotal));
	}

	//if(_totalAbundance == 1){
	//	cout << X2j << endl;
	//}
	if(X2j <= (_nbBanks-1)*1.5) return false;
	*/

	//cout << X2j << endl;


	//cout << kmer.toString(31) << endl;

	//cout << endl;

	//if(_progress){ //Simka_min
	//	_localStats->_nbSolidKmers += 1;
	//	computeStats(counts);
	//}
	//else{

		if(isAbundanceThreshold){

			if(!isSolidVector(counts))
				return false;

			/*
			cout << endl;
			for(size_t i=0; i<counts.size(); i++)
				cout << counts[i] << " ";
			cout << endl;*/
			//cout << _abundanceThreshold.first << " " << _abundanceThreshold.second << endl;

			for(size_t i=0; i<counts.size(); i++){

				if(counts[i] >= _abundanceThreshold.first && counts[i] <= _abundanceThreshold.second)
					_solidCounts[i] = counts[i];
				else
					_solidCounts[i] = 0;
			}


			//for(size_t i=0; i<counts.size(); i++)
			//	cout << _solidCounts[i] << " ";
			//cout << endl;

			computeStats(_solidCounts);

			//computeStats(counts);
		}
		else{
			computeStats(counts);

		}

		_localStats->_nbSolidKmers += 1;


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
			_localStats->_chord_N2[i] += pow(abundanceI, 2);
		}


		for(size_t j=i+1; j<counts.size(); j++){
			CountNumber abundanceJ = counts[j];

			/*
			if(_stats._distanceParams._computeBrayCurtis)
				_localStats->_brayCurtisNumerator[i][j] += abs(abundanceI - abundanceJ);

			if(_stats._distanceParams._computeChord)
				_localStats->_chord_NiNj[i][j] += abundanceI * abundanceJ;

			if(_stats._distanceParams._computeHellinger)
				_localStats->_hellinger_SqrtNiNj[i][j] += sqrt(abundanceI * abundanceJ);

			if(_stats._distanceParams._computeCanberra){
				if(abundanceI + abundanceJ > 0){
					_localStats->_canberra[i][j] += pow((abundanceI - abundanceJ) / (abundanceI + abundanceJ), 2);
				}
			}

			if(_stats._distanceParams._computeKulczynski)
				_localStats->_kulczynski_minNiNj[i][j] += min(abundanceI, abundanceJ);*/


			if(abundanceI + abundanceJ > 0){
				_localStats->_canberra[i][j] += abs(abundanceI - abundanceJ) / (abundanceI + abundanceJ);
				_localStats->_brayCurtisNumerator[i][j] += abs(abundanceI - abundanceJ);
			}

			if(abundanceI && abundanceJ){
				_localStats->_matrixNbSharedKmers[i][j] += abundanceI;
				_localStats->_matrixNbSharedKmers[j][i] += abundanceJ;
				_localStats->_matrixNbDistinctSharedKmers[i][j] += 1;

				_localStats->_chord_NiNj[i][j] += abundanceI * abundanceJ;
				_localStats->_hellinger_SqrtNiNj[i][j] += sqrt(abundanceI * abundanceJ);
				_localStats->_kulczynski_minNiNj[i][j] += min(abundanceI, abundanceJ);
			}

			/*
			if(abundanceI && abundanceJ){
				_localStats->_matrixNbSharedKmers[i][j] += abundanceI;
				_localStats->_matrixNbSharedKmers[j][i] += abundanceJ;
				_localStats->_matrixNbDistinctSharedKmers[i][j] += 1;
				//_localStats->_matrixNbDistinctSharedKmers[j][i] += 1;
				//updateKullbackLeibler(i, abundanceI, j, abundanceJ);
			}*/

		}

	}

	/*
	if(_stats._distanceParams._computeBrayCurtis){
		for(size_t i=0; i<counts.size(); i++){
			for(size_t j=i+1; j<counts.size(); j++){
				_localStats->_brayCurtisNumerator[i][j] += abs(counts[i] - counts[j]);
			}
		}
	}

	if(_stats._distanceParams._computeChord){
		for(size_t i=0; i<counts.size(); i++){
			_localStats->_chord_N2[i] += pow(counts[i], 2);
			for(size_t j=i+1; j<counts.size(); j++){
				_localStats->_chord_NiNj[i][j] += counts[i] * counts[j];
			}
		}
	}

	if(_stats._distanceParams._computeHellinger){
		for(size_t i=0; i<counts.size(); i++){
			for(size_t j=i+1; j<counts.size(); j++){
				_localStats->_hellinger_SqrtNiNj[i][j] += sqrt(counts[i] * counts[j]);
			}
		}
	}

	if(_stats._distanceParams._computeCanberra){
		for(size_t i=0; i<counts.size(); i++){
			CountNumber abundanceI = counts[i];
			for(size_t j=i+1; j<counts.size(); j++){
				CountNumber abundanceJ = counts[j];
				if(abundanceI + abundanceJ > 0)
					_localStats->_canberra[i][j] += abs(abundanceI - abundanceJ) / (abundanceI + abundanceJ);
			}
		}

	}

	if(_stats._distanceParams._computeKulczynski){
		for(size_t i=0; i<counts.size(); i++){
			for(size_t j=i+1; j<counts.size(); j++){
				_localStats->_kulczynski_minNiNj[i][j] += min(counts[i], counts[j]);
			}
		}
	}*/

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

/*
template<size_t span>
void SimkaCountProcessor<span>::updateBrayCurtis(int bank1, CountNumber abundance1, int bank2, CountNumber abundance2){


	//_localStats->_brayCurtisNumerator[bank1][bank2] += abs(abundance1-abundance2);
	_localStats->_brayCurtisNumerator[bank1][bank2] += min(abundance1, abundance2);
	_localStats->_brayCurtisNumerator[bank2][bank1] += min(abundance1, abundance2);
}*/


template<size_t span>
double SimkaCountProcessor<span>::getShannonIndex(const Type&  kmer){
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





























template<size_t span>
SimkaAlgorithm<span>::SimkaAlgorithm(IProperties* options)
:
Algorithm("simka", -1, options)
//_progress (0), _tmpPartitionsStorage(0), _tmpPartitions(0)
{


	_options = options;
	_stats = 0;
	//_simkaDistance = 0;
	_banks = 0;
	_processor = 0;




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


	//_banksInputFilename = _inputFilename + "_dsk_dataset_temp__";


}

template<size_t span>
SimkaAlgorithm<span>::~SimkaAlgorithm() {
}



template<size_t span>
void SimkaAlgorithm<span>::execute() {

	if(!setup()) return;
	if(!isInputValid()) return;

	createBank();

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
bool SimkaAlgorithm<span>::setup() {

	if(! createDirs() ) return false;

	try{
		layoutInputFilename();
	}
	catch (Exception& e){
		cout << "Syntax error in input file" << endl;
		return false;
	}

	_nbBanks = _bankNames.size();
	computeMaxReads();

	return true;
}


template<size_t span>
void SimkaAlgorithm<span>::parseArgs() {

	_maxMemory = _options->getInt(STR_MAX_MEMORY);
    _nbCores = _options->getInt(STR_NB_CORES);
	_inputFilename = _options->getStr(STR_URI_INPUT);
	_outputDir = _options->get(STR_URI_OUTPUT) ? _options->getStr(STR_URI_OUTPUT) : "./";
	_outputDirTemp = _options->get(STR_URI_OUTPUT_TMP) ? _options->getStr(STR_URI_OUTPUT_TMP) : "./";
	_kmerSize = _options->getInt(STR_KMER_SIZE);
	_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
	_abundanceThreshold.second = _options->getInt(STR_KMER_ABUNDANCE_MAX);
	_soliditySingle = _options->get(STR_SIMKA_SOLIDITY_PER_DATASET);
	//_nbMinimizers = _options->getInt(STR_KMER_PER_READ);
	//_maxDisk = getInput()->getInt(STR_MAX_DISK);

	//read filter
	_maxNbReads = _options->getInt(STR_SIMKA_MAX_READS);
	_minReadSize = _options->getInt(STR_SIMKA_MIN_READ_SIZE);
	_minReadShannonIndex = _options->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
	_minReadShannonIndex = std::max(_minReadShannonIndex, 0.0);
	_minReadShannonIndex = std::min(_minReadShannonIndex, 2.0);

	_minKmerShannonIndex = _options->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
	_minKmerShannonIndex = std::max(_minKmerShannonIndex, 0.0);
	_minKmerShannonIndex = std::min(_minKmerShannonIndex, 2.0);

}

template<size_t span>
bool SimkaAlgorithm<span>::createDirs(){

	if(!System::file().doesExist(_outputDir)){
		int ok = System::file().mkdir(_outputDir, -1);
		if(ok != 0){
	        std::cout << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
	        return false;
		}
	}

	_outputDirTemp = _outputDirTemp;

	if(!System::file().doesExist(_outputDirTemp)){
		int ok = System::file().mkdir(_outputDirTemp, -1);
		if(ok != 0){
	        std::cout << "Error: can't create output temp directory (" << _outputDirTemp << ")" << std::endl;
	        return false;
		}
	}

	_outputDirTemp = System::file().getRealPath(_outputDirTemp);
	_outputDirTemp += "/simka_output_temp/";
	System::file().mkdir(_outputDirTemp, -1);

	_options->setStr(STR_URI_OUTPUT_TMP, _outputDirTemp);
	System::file().mkdir(_outputDirTemp + "/input/", -1);

	return true;
}

template<size_t span>
void SimkaAlgorithm<span>::layoutInputFilename(){

	if(_options->getInt(STR_VERBOSE) != 0){
		cout << endl << "Creating input" << endl;
	}

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


	if(_options->getInt(STR_VERBOSE) != 0){
		cout << "\tNb input datasets: " << _bankNames.size() << endl;
		cout << endl;
	}

}


template<size_t span>
bool SimkaAlgorithm<span>::isInputValid(){

	string inputDir = _outputDirTemp + "/input/";

	for (size_t i=0; i<_nbBanks; i++){

		try{
			IBank* bank = Bank::open(inputDir + _bankNames[i]);
			LOCAL(bank);
		}
		catch (Exception& e){
			cout << "Can't open dataset: " << _bankNames[i] << endl;
			return false;
		}

	}

	return true;

}

template<size_t span>
void SimkaAlgorithm<span>::computeMaxReads(){

	string inputDir = _outputDirTemp + "/input/";

	//if(_maxNbReads != 0){
	//	return;
	//}

	if(_maxNbReads == 0){
		if(_options->getInt(STR_VERBOSE) != 0)
			cout << "-maxNbReads is not defined. Simka will estimating it..." << endl;
	}

	u_int64_t totalReads = 0;
	u_int64_t minReads = -1;
	u_int64_t maxReads = 0;
	for (size_t i=0; i<_nbBanks; i++){

		IBank* bank = Bank::open(inputDir + _bankNames[i]);
		LOCAL(bank);
		u_int64_t nbReads = bank->estimateNbItems();
		nbReads /= _nbBankPerDataset[i];
		totalReads += nbReads;
		if(nbReads < minReads){
			minReads = nbReads;
			//_smallerBankId = _bankNames[i];
		}
		if(nbReads > maxReads){
			maxReads = nbReads;
			_largerBankId = _bankNames[i];
		}

	}

	u_int64_t meanReads = totalReads / _nbBanks;

	if(_options->getInt(STR_VERBOSE) != 0){
		cout << "Smaller datasets contains: " << minReads << " reads" << endl;
		cout << "Larger datasets contains: " << maxReads << " reads" << endl;
		cout << "Whole datasets contains a mean of: " << meanReads << " reads" << endl;
	}

	if(_maxNbReads == 0){
		_maxNbReads = (minReads + meanReads) / 2;
		if(_options->getInt(STR_VERBOSE) != 0){
			cout << "Simka will use: " << _maxNbReads << " reads per dataset"<< endl << endl;
		}
	}
	else if(_maxNbReads == -1){
		if(_options->getInt(STR_VERBOSE) != 0)
			cout << "Simka will use all reads"<< endl << endl;
		_maxNbReads = 0;
	}

}

/*

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

			bankFileContents += subBankFilename + "\n";
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


}*/


template<size_t span>
void SimkaAlgorithm<span>::createBank(){

	IBank* bank = Bank::open(_banksInputFilename);

	SimkaSequenceFilter sequenceFilter(_minReadSize, _minReadShannonIndex);

	_banks = new SimkaBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, _nbBankPerDataset, _maxNbReads);

}

template<size_t span>
void SimkaAlgorithm<span>::count(){


	SimkaDistanceParam distanceParams(_options);

	_stats = new SimkaStatistics(_nbBanks, distanceParams);

	SortingCountAlgorithm<span> sortingCount (_banks, _options);

	// We create a custom count processor and give it to the sorting count algorithm
	_processor = new SimkaCountProcessor<span> (*_stats, _nbBanks, _kmerSize, _abundanceThreshold, _solidKind, _soliditySingle, _minKmerShannonIndex);
	_processor->use();
	sortingCount.addProcessor (_processor);

	// We launch the algorithm
	sortingCount.execute();


}









template<size_t span>
void SimkaAlgorithm<span>::outputMatrix(){
	_stats->outputMatrix(_outputDir, _bankNames);
}




template<size_t span>
void SimkaAlgorithm<span>::print(){

	cout << "Output folder:   " << _outputDir << endl;



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
