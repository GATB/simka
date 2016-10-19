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
SimkaAlgorithm<span>::SimkaAlgorithm(IProperties* options)
:
Algorithm("simka", -1, options)
//_progress (0), _tmpPartitionsStorage(0), _tmpPartitions(0)
{


	_options = options;
	_stats = 0;
	//_simkaDistance = 0;
	_banks = 0;
	//_processor = 0;




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

	/*
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

	clear();*/
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

	return true;
}


template<size_t span>
void SimkaAlgorithm<span>::parseArgs() {

	_computeSimpleDistances = _options->get(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
	_computeComplexDistances = _options->get(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);
	_keepTmpFiles = _options->get(STR_SIMKA_KEEP_TMP_FILES);
	_maxMemory = _options->getInt(STR_MAX_MEMORY);
    _nbCores = _options->getInt(STR_NB_CORES);
	_inputFilename = _options->getStr(STR_URI_INPUT);
	_outputDir = _options->get(STR_URI_OUTPUT) ? _options->getStr(STR_URI_OUTPUT) : "./";
	_outputDirTemp = _options->get(STR_URI_OUTPUT_TMP) ? _options->getStr(STR_URI_OUTPUT_TMP) : "./";
	_kmerSize = _options->getInt(STR_KMER_SIZE);
	_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
	_abundanceThreshold.second = min((u_int64_t)_options->getInt(STR_KMER_ABUNDANCE_MAX), (u_int64_t)(999999999));

	//cout << _options->getInt(STR_KMER_ABUNDANCE_MAX) << endl;
	//cout << _abundanceThreshold.second << endl;
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

	if(!System::file().doesExist(_inputFilename)){
		cerr << "ERROR: Input filename does not exist" << endl;
		exit(1);
	}
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
			cerr << "ERROR: Can't open dataset: " << _bankNames[i] << endl;
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
	u_int64_t meanReads = 0;

	if(_maxNbReads == 0 || _options->get(STR_SIMKA_COMPUTE_DATA_INFO)){

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

		meanReads = totalReads / _nbBanks;

		if(_options->getInt(STR_VERBOSE) != 0){
			cout << "Smaller sample contains: " << minReads << " reads" << endl;
			cout << "Larger sample contains: " << maxReads << " reads" << endl;
			cout << "Whole dataset contains a mean of: " << meanReads << " reads" << endl << endl;
		}
	}


	if(_maxNbReads == 0){
		_maxNbReads = (minReads + meanReads) / 2;
		if(_options->getInt(STR_VERBOSE) != 0){
			cout << "Reads per sample used up to: " << _maxNbReads << endl << endl;
		}
	}
	else if(_maxNbReads == -1){
		if(_options->getInt(STR_VERBOSE) != 0)
			cout << "Reads per sample used: all"<< endl << endl;
		_maxNbReads = 0;
	}
	else{
		if(_options->getInt(STR_VERBOSE) != 0){
			cout << "Reads per sample used up to: " << _maxNbReads << endl << endl;
		}
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

	/*
	//SimkaDistanceParam distanceParams(_options);

	_stats = new SimkaStatistics(_nbBanks, _computeSimpleDistances, _computeComplexDistances, _outputDirTemp, _bankNames);

	SortingCountAlgorithm<span> sortingCount (_banks, _options);

	// We create a custom count processor and give it to the sorting count algorithm
	vector<u_int64_t> dummyVec;
	_processor = new SimkaCountProcessor<span> (*_stats, _nbBanks, _kmerSize, _abundanceThreshold, _solidKind, _soliditySingle, _minKmerShannonIndex);
	_processor->use();
	sortingCount.addProcessor (_processor);

	// We launch the algorithm
	sortingCount.execute();
	*/

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
    //if(_processor) _processor->forget();

    for(size_t i=0; i<_tempFilenamesToDelete.size(); i++){
    	System::file().remove(_tempFilenamesToDelete[i]);
    }

    if(_stats) delete _stats;
    //if(_simkaDistance) delete _simkaDistance;
	//_banks->remove();
	//delete _processor;
}
