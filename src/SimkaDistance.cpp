/*
 * SimkaDistance.cpp
 *
 *  Created on: 17 juin 2015
 *      Author: gbenoit
 */

#include "SimkaDistance.hpp"





SimkaStatistics::SimkaStatistics(size_t nbBanks){

	_nbBanks = nbBanks;

	_nbKmers = 0;
	_nbDistinctKmers = 0;
	_nbSolidKmers = 0;
	_nbErroneousKmers = 0;

	//_abundanceMin = abundanceMin;
	//_mutex = mutex;
	//_outputDir = outputDir;

	_nbSolidDistinctKmersPerBank.resize(_nbBanks, 0);
	_nbSolidKmersPerBank.resize(_nbBanks, 0);
	_nbKmersPerBank.resize(_nbBanks, 0);

	_nbDistinctKmersSharedByBanksThreshold.resize(_nbBanks, 0);
	_nbKmersSharedByBanksThreshold.resize(_nbBanks, 0);

	_matrixNbDistinctSharedKmers.resize(_nbBanks);
	_matrixNbSharedKmers.resize(_nbBanks);
	_brayCurtisNumerator.resize(_nbBanks);
	_kullbackLeibler.resize(_nbBanks);
	for(int i=0; i<_nbBanks; i++){
		_matrixNbDistinctSharedKmers[i].resize(nbBanks, 0);
		_matrixNbSharedKmers[i].resize(nbBanks, 0);
		_brayCurtisNumerator[i].resize(nbBanks, 0);
		_kullbackLeibler[i].resize(nbBanks, 0);
	}

}


SimkaStatistics& SimkaStatistics::operator+=  (const SimkaStatistics& other){
	_nbKmers += other._nbKmers;
	_nbDistinctKmers += other._nbDistinctKmers;
	_nbSolidKmers += other._nbSolidKmers;
	_nbErroneousKmers += other._nbErroneousKmers;

	for(size_t i=0; i<_nbBanks; i++){
		_nbKmersPerBank[i] += other._nbKmersPerBank[i];
		_nbSolidDistinctKmersPerBank[i] += other._nbSolidDistinctKmersPerBank[i];
		_nbSolidKmersPerBank[i] += other._nbSolidKmersPerBank[i];
		_nbDistinctKmersSharedByBanksThreshold[i] += other._nbDistinctKmersSharedByBanksThreshold[i];
		_nbKmersSharedByBanksThreshold[i] += other._nbKmersSharedByBanksThreshold[i];
	}


	for(size_t i=0; i<_nbBanks; i++){
		for(size_t j=0; j<_nbBanks; j++){
			_matrixNbDistinctSharedKmers[i][j] += other._matrixNbDistinctSharedKmers[i][j];
			_matrixNbSharedKmers[i][j] += other._matrixNbSharedKmers[i][j];
			_brayCurtisNumerator[i][j] += other._brayCurtisNumerator[i][j];
			_kullbackLeibler[i][j] += other._kullbackLeibler[i][j];
		}
	}

	return *this;
}

void SimkaStatistics::print(){

	//cout.precision(4);
    cout << endl << endl;

    //return;

    u_int64_t solidAbundance = 0;
    //for(int i=0; i<_nbSolidKmersPerBankAbundance.size(); i++)
    //	solidAbundance += _nbSolidKmersPerBankAbundance[i];
    for(int i=0; i<_nbKmersSharedByBanksThreshold.size(); i++)
    	solidAbundance += _nbKmersSharedByBanksThreshold[i];

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
	    cout << "\t\t\tDistinct:    " << _nbDistinctKmersSharedByBanksThreshold[i] << "    ";
	    if(_nbSolidKmers > 0){
		    cout << (_nbDistinctKmersSharedByBanksThreshold[i]*100) / (float)_nbSolidKmers << "%";
	    }
	    else{
	    	cout << "0%";
	    }

	    cout << endl;
	    cout << "\t\t\tAbundance:    " << _nbKmersSharedByBanksThreshold[i] << "    ";
	    if(solidAbundance > 0){
	    	cout << (_nbKmersSharedByBanksThreshold[i]*100) / (float)solidAbundance << "%";
	    }
	    else{
	    	cout << "0%";
	    }
	    if(_nbDistinctKmersSharedByBanksThreshold[i] > 0){
	    	cout << endl;
		    cout << "\t\t\tMean abundance per bank: " << _nbKmersSharedByBanksThreshold[i] / _nbDistinctKmersSharedByBanksThreshold[i] / (float) _nbBanks;
	    }

	    cout << endl;
    }

    //cout << endl;
    //cout << "Nb kmers in all banks (max/min > 10): " << _nbKmersInCoupleBankSupRatio << "    " << (_nbKmersInCoupleBankSupRatio*100) / (float)_nbSolidKmers << "%" <<  endl;


   cout << endl << endl;
}










SimkaDistance::SimkaDistance(SimkaStatistics& stats) : _stats(stats){
	_nbBanks = _stats._nbBanks;
}

//SimkaDistance::~SimkaDistance() {
	// TODO Auto-generated destructor stub
//}

/*

void SimkaDistance::outputMatrix(){


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
	    	matrixNormalized[i][j] = (100.0 * (_processor->_matrixNbDistinctSharedKmers[i][j] + _processor->_matrixNbDistinctSharedKmers[j][i])) / (_processor->_nbSolidDistinctKmersPerBank[i] + _processor->_nbSolidDistinctKmersPerBank[j]);
	    	matrixPercentage[i][j] = (100.0 * (_processor->_matrixNbDistinctSharedKmers[i][j])) / (_processor->_nbSolidDistinctKmersPerBank[i]);
	    	matrixAbundanceNormalized[i][j] = (100.0 * (_processor->_matrixNbSharedKmers[i][j] + _processor->_matrixNbSharedKmers[j][i])) / (_processor->_nbSolidKmersPerBank[i] + _processor->_nbSolidKmersPerBank[j]);
	    	matrixAbundancePercentage[i][j] = (100.0 * (_processor->_matrixNbSharedKmers[i][j])) / (_processor->_nbSolidKmersPerBank[i]);
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
	_matDksPercFilename = "mat_dks_asym" + strKmerSize + strAbMin + strAbMax + ".csv";
	_matAksNormFilename = "mat_aks_norm" + strKmerSize + strAbMin + strAbMax + ".csv";
	_matAksPercFilename = "mat_aks_asym" + strKmerSize + strAbMin + strAbMax + ".csv";

    dumpMatrix(_matDksNormFilename, matrixNormalized);
    dumpMatrix(_matDksPercFilename, matrixPercentage);
    dumpMatrix(_matAksNormFilename, matrixAbundanceNormalized);
    dumpMatrix(_matAksPercFilename, matrixAbundancePercentage);


}
*/

vector<vector<float> > SimkaDistance::createSquaredMatrix(int n){
    vector<vector<float> > matrix;

    matrix.resize(n);
    for(int i=0; i<n; i++)
    	matrix[i].resize(n, 0);

    return matrix;

}

//Presence/absence Sorensen 2c/(S1+S2)
vector<vector<float> > SimkaDistance::getMatrixSorensen(SIMKA_MATRIX_TYPE type){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

    if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
    	for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (_stats._matrixNbDistinctSharedKmers[i][j])) / (_stats._nbSolidDistinctKmersPerBank[i]);
    }
    else if(type == SIMKA_MATRIX_TYPE::NORMALIZED){
        for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (2*_stats._matrixNbDistinctSharedKmers[i][j])) / (_stats._nbSolidDistinctKmersPerBank[i] + _stats._nbSolidDistinctKmersPerBank[j]);
    }

    return matrix;
}

//Presence/absence Jaccard |AnB| / |AuB|
vector<vector<float> > SimkaDistance::getMatrixJaccard(){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

    //if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
    //    for(int i=0; i<_nbBanks; i++)
    //	    for(int j=0; j<_nbBanks; j++)
    //	    	matrix[i][j] = (100.0 * (_stats._matrixNbDistinctSharedKmers[i][j])) / (_stats._nbSolidDistinctKmersPerBank[i]);
    //}
    //else if(type == SIMKA_MATRIX_TYPE::NORMALIZED){

        for(int i=0; i<_nbBanks; i++){
    	    for(int j=0; j<_nbBanks; j++){
    	    	u_int64_t intersection_ = _stats._matrixNbDistinctSharedKmers[i][j];
    	    	u_int64_t union_ = _stats._nbSolidDistinctKmersPerBank[i] + _stats._nbSolidDistinctKmersPerBank[j] - _stats._matrixNbDistinctSharedKmers[i][j];
    	    	matrix[i][j] = (100.0 * intersection_) / union_;
    	    }
		}

    return matrix;
}

//abundance
vector<vector<float> > SimkaDistance::getMatrixAKS(SIMKA_MATRIX_TYPE type){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

    if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
        for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (_stats._matrixNbSharedKmers[i][j])) / (_stats._nbSolidKmersPerBank[i]);
    }
    else if(type == SIMKA_MATRIX_TYPE::NORMALIZED){

        for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (_stats._matrixNbSharedKmers[i][j] + _stats._matrixNbSharedKmers[j][i])) / (_stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j]);
    }

    return matrix;
}

//Abundance: bray curtis
vector<vector<float> > SimkaDistance::getMatrixBrayCurtis(){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

	for(int i=0; i<_nbBanks; i++)
		for(int j=0; j<_nbBanks; j++)
			matrix[i][j] = (100.0 * (2*_stats._brayCurtisNumerator[i][j])) / (_stats._nbSolidKmersPerBank[i]+_stats._nbSolidKmersPerBank[j]);

    return matrix;
}

//Abundance: Kullback Leibler
vector<vector<float> > SimkaDistance::getMatrixKullbackLeibler(){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

	for(int i=0; i<_nbBanks; i++)
		for(int j=0; j<_nbBanks; j++)
			matrix[i][j] = _stats._kullbackLeibler[i][j];// / (_stats._nbSolidKmersPerBank[i]); // (100.0 * (2*_stats._kullbackLeibler[i][j])) / (_stats._nbSolidKmersPerBank[i]+_stats._nbSolidKmersPerBank[j]);

    return matrix;
}
