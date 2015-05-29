/*
 * SimkaAlgorithm.cpp
 *
 *  Created on: 18 mai 2015
 *      Author: gbenoit
 */

#include "SimkaAlgorithm.hpp"
























template<size_t span>
SimkaCountProcessor<span>::SimkaCountProcessor (size_t nbBanks, const pair<size_t, size_t>& abundanceThreshold){
	// We configure the vector for the N.(N+1)/2 possible pairs
	//_countTotal.resize (_nbBanks*(_nbBanks+1)/2);

	_nbBanks = nbBanks;
	_abundanceThreshold = abundanceThreshold;

	_nbKmers = 0;
	_nbDistinctKmers = 0;
	_nbSolidKmers = 0;
	_nbErroneousKmers = 0;

	//_abundanceMin = abundanceMin;
	//_mutex = mutex;
	//_outputDir = outputDir;

	_nbSolidKmersPerBank.resize(_nbBanks, 0);
	_nbSolidKmersPerBankAbundance.resize(_nbBanks, 0);
	_nbKmersPerBank.resize(_nbBanks, 0);

	_nbKmersSharedByBanksThreshold.resize(_nbBanks, 0);
	_nbKmersAbundanceSharedByBanksThreshold.resize(_nbBanks, 0);

	_matrixSharedKmers.resize(_nbBanks);
	_matrixSharedAbundanceKmers.resize(_nbBanks);
	for(int i=0; i<_nbBanks; i++){
		_matrixSharedKmers[i].resize(nbBanks, 0);
		_matrixSharedAbundanceKmers[i].resize(nbBanks, 0);
	}

}

template<size_t span>
SimkaCountProcessor<span>::~SimkaCountProcessor () {

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

	_nbKmers += clone->_nbKmers;
	_nbDistinctKmers += clone->_nbDistinctKmers;
	_nbSolidKmers += clone->_nbSolidKmers;
	_nbErroneousKmers += clone->_nbErroneousKmers;

	for(size_t i=0; i<_nbBanks; i++){
		_nbKmersPerBank[i] += clone->_nbKmersPerBank[i];
		_nbSolidKmersPerBank[i] += clone->_nbSolidKmersPerBank[i];
		_nbSolidKmersPerBankAbundance[i] += clone->_nbSolidKmersPerBankAbundance[i];
		_nbKmersSharedByBanksThreshold[i] += clone->_nbKmersSharedByBanksThreshold[i];
		_nbKmersAbundanceSharedByBanksThreshold[i] += clone->_nbKmersAbundanceSharedByBanksThreshold[i];
	}


	for(size_t i=0; i<_nbBanks; i++){
		for(size_t j=0; j<_nbBanks; j++){
			_matrixSharedKmers[i][j] += clone->_matrixSharedKmers[i][j];
			_matrixSharedAbundanceKmers[i][j] += clone->_matrixSharedAbundanceKmers[i][j];
		}
	}

}

template<size_t span>
bool SimkaCountProcessor<span>::isSolid(const CountVector& counts){

	bool isSolid = false;

	for(size_t i=0; i<counts.size(); i++){

		CountNumber abundance = counts[i];
		if(abundance == 0) continue;

		if(abundance >= _abundanceThreshold.first && abundance <= _abundanceThreshold.second){
			return true;
		}

		//nbBanks += 1;
		//if(nbBanks > 1){
		//	isSolid = true;
		//	break;
		//}
	}

	return false;

}

template<size_t span>
bool SimkaCountProcessor<span>::process (size_t partId, const Type& kmer, const CountVector& counts, CountNumber sum){

	_nbDistinctKmers += 1;

	for(size_t i=0; i<counts.size(); i++){
		_nbKmers += counts[i];
		_nbKmersPerBank[i] += counts[i];
	}


	if(isSolid(counts))
		_nbSolidKmers += 1;
	else
		return false;


	int nbBanksThatHaveKmer = 0;
	//vector<bool> hasBankKmer(counts.size(), false);
	u_int64_t totalAbundance = 0;

	//pair<u_int32_t, u_int32_t> kmerInBankCoupleAbundance(-1, -1);

	for(size_t i=0; i<counts.size(); i++){

		CountNumber abundance = counts[i];
		//if(abundance < _abundanceMin) continue;

		if(abundance != 0){
			totalAbundance += abundance;
			nbBanksThatHaveKmer += 1;
			_nbSolidKmersPerBank[i] += 1;
			_nbSolidKmersPerBankAbundance[i] += abundance;
			//hasBankKmer[i] = true;

			//if(kmerInBankCoupleAbundance.first == -1)
			//	kmerInBankCoupleAbundance.first = abundance;
			//else if(kmerInBankCoupleAbundance.second == -1)
			//	kmerInBankCoupleAbundance.second = abundance;
		}


	}


	for(size_t i=0; i<counts.size(); i++){
		for(size_t j=0; j<counts.size(); j++){
			if(counts[i] && counts[j]){
				_matrixSharedAbundanceKmers[i][j] += counts[i];
				_matrixSharedKmers[i][j] += 1;
			}

		}
	}

	_nbKmersSharedByBanksThreshold[nbBanksThatHaveKmer-1] += 1;
	_nbKmersAbundanceSharedByBanksThreshold[nbBanksThatHaveKmer-1] += totalAbundance;

	if(totalAbundance == 1){
		//if( == 1){
			_nbErroneousKmers += 1;
		//}
	}
	//else if(nbBanksThatHaveKmer == counter.size()){
	//}


	return true;
}

template<size_t span>
void SimkaCountProcessor<span>::print(){

	cout.precision(4);
    cout << endl << endl;

    //return;

    u_int64_t solidAbundance = 0;
    for(int i=0; i<_nbKmersAbundanceSharedByBanksThreshold.size(); i++)
    	solidAbundance += _nbKmersAbundanceSharedByBanksThreshold[i];

    cout << "Statistics on kmer intersections:" << endl;
    cout << "\tNb kmers: " << _nbKmers << "    " << _nbKmers / 1000000 << " M" << "    " << _nbKmers / 1000000000 << " G" << endl;
    cout << endl;

    cout << "\tNb distinct kmers: " << _nbDistinctKmers << "    " << _nbDistinctKmers / 1000000 << " M" << "    " << _nbDistinctKmers / 1000000000 << " G" << "    " << (100*_nbDistinctKmers)/(float)_nbKmers << "%" << endl;
    cout << "\tNb solid kmers: " << _nbSolidKmers << "    " << _nbSolidKmers / 1000000 << " M" << "    " << _nbSolidKmers / 1000000000 << " G" << "    " << (100*_nbSolidKmers)/(float)_nbDistinctKmers << "% distinct" << "       " << (100*solidAbundance) / (double)_nbKmers << "% abundance" << endl;
    //for(int i=0; i<_nbBanks; i++){
	    //cout << "Nb kmers (M) " << i <<  ": " << _nbSolidKmersPerBank[i] << endl << endl;
    //}

    cout << endl;
    cout << "\tPotentially erroneous (Kmers appearing only one time in a single bank): " << endl;
    cout << "\t\t" << _nbErroneousKmers << "    " << _nbErroneousKmers / 1000000 << " M" << "    " << _nbErroneousKmers / 1000000000 << " G" << "    " << (100*_nbErroneousKmers)/(float)_nbDistinctKmers << "% distinct" << "      " << (100*_nbErroneousKmers)/(float)_nbKmers << "% abundance" << endl;

    cout << endl;
    cout << "\tKmer shared by T banks :" << endl;

    for(int i=0; i<_nbBanks; i++){
	    cout << "\t\tShared by " << i+1 <<  " banks:";

	    cout << endl;
	    cout << "\t\t\tDistinct:    " << _nbKmersSharedByBanksThreshold[i] << "    ";
	    if(_nbSolidKmers > 0){
		    cout << (_nbKmersSharedByBanksThreshold[i]*100) / (float)_nbSolidKmers << "%";
	    }
	    else{
	    	cout << "0%";
	    }

	    cout << endl;
	    cout << "\t\t\tAbundance:    " << _nbKmersAbundanceSharedByBanksThreshold[i] << "    ";
	    if(solidAbundance > 0){
	    	cout << (_nbKmersAbundanceSharedByBanksThreshold[i]*100) / (float)solidAbundance << "%";
	    }
	    else{
	    	cout << "0%";
	    }
	    if(_nbKmersSharedByBanksThreshold[i] > 0){
	    	cout << endl;
		    cout << "\t\t\tMean abundance per bank: " << _nbKmersAbundanceSharedByBanksThreshold[i] / _nbKmersSharedByBanksThreshold[i] / (float) _nbBanks;
	    }

	    cout << endl;
    }

    //cout << endl;
    //cout << "Nb kmers in all banks (max/min > 10): " << _nbKmersInCoupleBankSupRatio << "    " << (_nbKmersInCoupleBankSupRatio*100) / (float)_nbSolidKmers << "%" <<  endl;


   cout << endl << endl;
}


































template<size_t span>
SimkaAlgorithm<span>::SimkaAlgorithm(IProperties* options) {

	_options = options;


	_inputFilename = _options->getStr(STR_URI_INPUT);
	_outputDir = _options->get(STR_URI_OUTPUT) ? _options->getStr(STR_URI_OUTPUT) : "./";
	_kmerSize = _options->getInt(STR_KMER_SIZE);
	_abundanceThreshold.first = _options->getInt(STR_KMER_ABUNDANCE_MIN);
	_abundanceThreshold.second =  _options->getInt(STR_KMER_ABUNDANCE_MAX);

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


	_banksInputFilename = _inputFilename + ".dsk_banks.____temp";


}

template<size_t span>
SimkaAlgorithm<span>::~SimkaAlgorithm() {
}


template<size_t span>
void SimkaAlgorithm<span>::execute() {

	layoutInputFilename();
	_banks = Bank::open(_banksInputFilename);
	count();

	if(_options->getInt(STR_VERBOSE) > 0)
	    _processor->print();

	outputMatrix();
	outputHeatmap();
	clear();

}

template<size_t span>
void SimkaAlgorithm<span>::layoutInputFilename(){

	IFile* inputFile = System::file().newFile(_inputFilename, "rb");
	IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

	inputFile->seeko(0, SEEK_END);
	u_int64_t size = inputFile->tell();
	inputFile->seeko(0, SEEK_SET);
	char buffer2[size];
	inputFile->fread(buffer2, size, size);
	string fileContents(buffer2, size);

	string line;
	string linePart;
	vector<string> linePartList;
	stringstream fileContentsStream(fileContents);

	string bankFileContents = "";

	while(getline(fileContentsStream, line)){

		stringstream lineStream(line);
		linePartList.clear();

		//stringstream test(fileContents);

		while(getline(lineStream, linePart, ' ')){

			if(linePart != ""){
				linePartList.push_back(linePart);
			}
		}

		string bankId = linePartList[0];
		string bankFilename = linePartList[1];

		_bankNames.push_back(bankId);
		bankFileContents += bankFilename + "\n";

		//cout << bankId << "            "  << bankFilename << endl;
		//string nbReadsOrCoverage = linePartList[2];

		//map<string, string> genomeInfos;
		//genomeInfos["filename"] = filename;
		//genomeInfos["nbReadsOrCoverage"] = nbReadsOrCoverage;

		//_genomeInfos[genomeId] = genomeInfos;
		//cout << genomeId << endl;

	}

	bankFileContents.pop_back(); // "remove last /n

	bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);

	delete inputFile;
	bankFile->flush();
	delete bankFile;

}


template<size_t span>
void SimkaAlgorithm<span>::count(){

	_nbBanks = _banks->getCompositionNb();

	SortingCountAlgorithm<span> sortingCount (_banks, _options);

	// We create a custom count processor and give it to the sorting count algorithm
	_processor = new SimkaCountProcessor<span> (_nbBanks, _abundanceThreshold);
	_processor->use();
	sortingCount.setProcessor (_processor);

	// We launch the algorithm
	sortingCount.execute();


}


template<size_t span>
void SimkaAlgorithm<span>::outputMatrix(){



    vector<vector<float> > matrixNormalized;
    vector<vector<float> > matrixPercentage;
    vector<vector<float> > matrixAbundanceNormalized;
    vector<vector<float> > matrixAbundancePercentage;

    matrixNormalized.resize(_nbBanks);
    matrixPercentage.resize(_nbBanks);
    matrixAbundanceNormalized.resize(_nbBanks);
    matrixAbundancePercentage.resize(_nbBanks);
    for(int i=0; i<_nbBanks; i++){
    	matrixNormalized[i].resize(_nbBanks, 0);
    	matrixPercentage[i].resize(_nbBanks, 0);
    	matrixAbundanceNormalized[i].resize(_nbBanks, 0);
    	matrixAbundancePercentage[i].resize(_nbBanks, 0);
    }


    for(int i=0; i<_nbBanks; i++){
	    for(int j=0; j<_nbBanks; j++){
	    	matrixNormalized[i][j] = (100.0 * (_processor->_matrixSharedKmers[i][j] + _processor->_matrixSharedKmers[j][i])) / (_processor->_nbSolidKmersPerBank[i] + _processor->_nbSolidKmersPerBank[j]);
	    	matrixPercentage[i][j] = (100.0 * (_processor->_matrixSharedKmers[i][j])) / (_processor->_nbSolidKmersPerBank[i]);
	    	matrixAbundanceNormalized[i][j] = (100.0 * (_processor->_matrixSharedAbundanceKmers[i][j] + _processor->_matrixSharedAbundanceKmers[j][i])) / (_processor->_nbSolidKmersPerBankAbundance[i] + _processor->_nbSolidKmersPerBankAbundance[j]);
	    	matrixAbundancePercentage[i][j] = (100.0 * (_processor->_matrixSharedAbundanceKmers[i][j])) / (_processor->_nbSolidKmersPerBankAbundance[i]);
	    }
    }


	char buffer[200];

	string strKmerSize = "_k";
	snprintf(buffer,200,"%llu",_kmerSize);
	strKmerSize += string(buffer);

    string strAbMax = "";
    if(_abundanceThreshold.second < 1000000){
    	snprintf(buffer,200,"%llu",_abundanceThreshold.second);
    	strAbMax += "_max" + string(buffer);
    }

	string strAbMin = "_min";
	snprintf(buffer,200,"%llu",_abundanceThreshold.first);
	strAbMin += string(buffer);

	_matDksNormFilename = "mat_dks_norm" + strKmerSize + strAbMin + strAbMax + ".csv";
	_matDksPercFilename = "mat_dks_perc" + strKmerSize + strAbMin + strAbMax + ".csv";
	_matAksNormFilename = "mat_aks_norm" + strKmerSize + strAbMin + strAbMax + ".csv";
	_matAksPercFilename = "mat_aks_perc" + strKmerSize + strAbMin + strAbMax + ".csv";

    dumpMatrix(_matDksNormFilename, matrixNormalized);
    dumpMatrix(_matDksPercFilename, matrixPercentage);
    dumpMatrix(_matAksNormFilename, matrixAbundanceNormalized);
    dumpMatrix(_matAksPercFilename, matrixAbundancePercentage);

}

template<size_t span>
void SimkaAlgorithm<span>::dumpMatrix(const string& outputFilename, vector<vector<float> >& matrix){

	char buffer[200];
	string str;

	for(int i=0; i<_nbBanks; i++){
		str += ";" + _bankNames[i];
		//str += ";" + datasetInfos[i]._name;
	}
	str += '\n';

	for(int i=0; i<matrix.size(); i++){

		str += _bankNames[i] + ";";
		//str += datasetInfos[i]._name + ";";
		for(int j=0; j<matrix.size(); j++){

			//snprintf(buffer,200,"%.2f", matrix[i][j]);
			snprintf(buffer,200,"%f", matrix[i][j]);
			str += string(buffer) + ";";

			//str += to_string(matrix[i][j]) + ";";
		}

		//matrixNormalizedStr.erase(matrixNormalizedStr.end()-1);
		str.pop_back(); //remove ; at the end of the line
		str += '\n';
	}


	gatb::core::system::IFile* file = gatb::core::system::impl::System::file().newFile(_outputDir + "/" + outputFilename, "wb");
	file->fwrite(str.c_str(), str.size(), 1);
	file->flush();
	delete file;

}

template<size_t span>
void SimkaAlgorithm<span>::outputHeatmap(){
	__outputHeatmap(_matDksPercFilename, _matDksNormFilename);
	__outputHeatmap(_matAksPercFilename, _matAksNormFilename);
}


template<size_t span>
void SimkaAlgorithm<span>::__outputHeatmap(const string& matrixPercFilename, const string& matrixNormFilename){

	string filename = matrixPercFilename.substr(0, matrixPercFilename.size()-4); //remove mat extension .csv
	vector<string> linePartList;
	string outputFilename = "heatmap";
	string part;
	//string linePart;
	//vector<string> linePartList;
	stringstream stream(filename);

	while(getline(stream, part, '_')){
		//cout << part << endl;
		linePartList.push_back(part);
	}
	linePartList.erase(linePartList.begin() + 0); //Remove 'mat' prefix
	linePartList.erase(linePartList.begin() + 1); //Remove 'norm' or 'perc'

	for(size_t i=0; i<linePartList.size(); i++){
		outputFilename += "_" + linePartList[i];
	}
	outputFilename += ".png";

	string command = "Rscript ./Rscripts/heatmap.r " + _outputDir + "/" + matrixPercFilename + " " + _outputDir + "/" + matrixPercFilename + " " + _outputDir + "/" + outputFilename;
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

}

template<size_t span>
void SimkaAlgorithm<span>::clear(){

	System::file().remove(_banksInputFilename);
    _processor->forget();
	//_banks->remove();
	//delete _processor;
}



















template class SimkaAlgorithm <KSIZE_1>;
template class SimkaAlgorithm <KSIZE_2>;
template class SimkaAlgorithm <KSIZE_3>;
template class SimkaAlgorithm <KSIZE_4>;
//template class SimkaAlgorithm<KSIZE_1>;




