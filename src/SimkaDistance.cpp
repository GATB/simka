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
	//_kullbackLeibler.resize(_nbBanks);
	for(int i=0; i<_nbBanks; i++){
		_matrixNbDistinctSharedKmers[i].resize(nbBanks, 0);
		_matrixNbSharedKmers[i].resize(nbBanks, 0);
		_brayCurtisNumerator[i].resize(nbBanks, 0);
		//_kullbackLeibler[i].resize(nbBanks, 0);
	}

	_speciesAbundancePerDataset.resize(_nbBanks);
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
			//_kullbackLeibler[i][j] += other._kullbackLeibler[i][j];
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

	u_int64_t a;
	u_int64_t b;
	u_int64_t c;

    _matrixBrayCurtis = createSquaredMatrix(_nbBanks);
    _matrixSymJaccardPresenceAbsence = createSquaredMatrix(_nbBanks);
    _matrixAsymJaccardPresenceAbsence = createSquaredMatrix(_nbBanks);
    _matrixSymJaccardAbundance = createSquaredMatrix(_nbBanks);
    _matrixAsymJaccardAbundance = createSquaredMatrix(_nbBanks);

    _matrixSymSorensen = createSquaredMatrix(_nbBanks);
    _matrixAsymSorensen = createSquaredMatrix(_nbBanks);
    //_matrixKullbackLeibler = createSquaredMatrix(_nbBanks);

	for(size_t i=0; i<_nbBanks; i++){
		//SpeciesAbundanceVectorType& X_i = _stats._speciesAbundancePerDataset[i];

		for(size_t j=i+1; j<_nbBanks; j++){
			//SpeciesAbundanceVectorType& X_j = _stats._speciesAbundancePerDataset[j];

			get_abc(i, j, a, b ,c);

			//Jaccard
			double jaccardSym = jaccardSimilarity(i, j, SYMETRICAL, PRESENCE_ABSENCE);
			_matrixSymJaccardPresenceAbsence[i][j] = jaccardSym;
			_matrixSymJaccardPresenceAbsence[j][i] = jaccardSym;

			_matrixAsymJaccardPresenceAbsence[i][j] = jaccardSimilarity(i, j, ASYMETRICAL, PRESENCE_ABSENCE);
			_matrixAsymJaccardPresenceAbsence[j][i] = jaccardSimilarity(j, i, ASYMETRICAL, PRESENCE_ABSENCE);

			jaccardSym = jaccardSimilarity(i, j, SYMETRICAL, ABUNDANCE);
			_matrixSymJaccardAbundance[i][j] = jaccardSym;
			_matrixSymJaccardAbundance[j][i] = jaccardSym;

			_matrixAsymJaccardAbundance[i][j] = jaccardSimilarity(i, j, ASYMETRICAL, ABUNDANCE);
			_matrixAsymJaccardAbundance[j][i] = jaccardSimilarity(j, i, ASYMETRICAL, ABUNDANCE);

			//Sorensen
			double sorensenSym = sorensenSimilarity(a, b, c, i, j, SYMETRICAL);
			_matrixSymSorensen[i][j] = sorensenSym;
			_matrixSymSorensen[j][i] = sorensenSym;

			_matrixAsymSorensen[i][j] = sorensenSimilarity(a, b, c, i, j, ASYMETRICAL);
			_matrixAsymSorensen[j][i] = sorensenSimilarity(a, b, c, j, i, ASYMETRICAL);



	        double bray = brayCurtisSimilarity(i,j);
			_matrixBrayCurtis[i][j] = bray;
			_matrixBrayCurtis[j][i] = bray;


			//double kl_val = kullbackLeibler(X_i,X_j);
			//_matrixKullbackLeibler[i][j] = kl_val;
			//_matrixKullbackLeibler[j][i] = kl_val;
		}

		/*
		for(size_t j=0; j<_nbBanks; j++){

			//Jaccard
			double jaccardAsym = jaccardSimilarity(i, j, ASYMETRICAL);
			_matrixSymJaccard[i][j] = jaccardAsym;

		}*/

	}
	/*
	    X_i = A[:,i]
	    for j in range(i+1,p):
	        X_j = A[:,j]
	        # getting a,b, and c
	        abc = get_abc(X_i,X_j)
	        # Jaccard distances
	        J = jaccard(abc)
	        Jac[i,j] = J
	        Jac[j,i] = J
	        # soerensen distances
	        S = soerensen(abc)
	        Sor[i,j] = S
	        Sor[j,i] = S
	        # Bray-Curtis distances
	        B = bc(X_i,X_j)
	        BC[i,j] = B
	        BC[j,i] = B
	        # Kullback-Leibler divergences
	        kl_val = kl(X_i,X_j)
	        KL[i,j] = kl_val
	        KL[j,i] = kl_val*/
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
    	matrix[i].resize(n, 100);

    return matrix;

}

/*
//Presence/absence Sorensen 2c/(S1+S2)
vector<vector<float> > SimkaDistance::getMatrixSorensen(SIMKA_MATRIX_TYPE type){

	return _matrixSorensen;

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
}*/

/*
//Presence/absence Jaccard |AnB| / |AuB|
vector<vector<float> > SimkaDistance::getMatrixJaccard(){

	return _matrixSymJaccard;
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
}*/

/*
//abundance
vector<vector<float> > SimkaDistance::getMatrixAKS(SIMKA_MATRIX_TYPE type){

    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

    if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
        for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (_stats._matrixNbSharedKmers[i][j])) / (_stats._nbSolidKmersPerBank[i]);
    }
    else if(type == SIMKA_MATRIX_TYPE::SYMETRICAL){

        for(int i=0; i<_nbBanks; i++)
    	    for(int j=0; j<_nbBanks; j++)
    	    	matrix[i][j] = (100.0 * (_stats._matrixNbSharedKmers[i][j] + _stats._matrixNbSharedKmers[j][i])) / (_stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j]);
    }

    return matrix;
}*/

/*
//Abundance: bray curtis
vector<vector<float> > SimkaDistance::getMatrixBrayCurtis(){

	return _matrixBrayCurtis;
    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

	for(int i=0; i<_nbBanks; i++)
		for(int j=0; j<_nbBanks; j++)
			//matrix[i][j] = (100.0 * (1 - ((_stats._brayCurtisNumerator[i][j]) / (float)(_stats._nbSolidKmersPerBank[i]+_stats._nbSolidKmersPerBank[j]))));
			matrix[i][j] = (100.0 * (2*_stats._brayCurtisNumerator[i][j])) / (_stats._nbSolidKmersPerBank[i]+_stats._nbSolidKmersPerBank[j]);

    return matrix;
}*/

//Abundance: Kullback Leibler
//vector<vector<float> > SimkaDistance::getMatrixKullbackLeibler(){

	//return _matrixKullbackLeibler;

	/*
    vector<vector<float> > matrix = createSquaredMatrix(_nbBanks);

	for(int i=0; i<_nbBanks; i++)
		for(int j=0; j<_nbBanks; j++)
			matrix[i][j] = _stats._kullbackLeibler[i][j];// / (_stats._nbSolidKmersPerBank[i]); // (100.0 * (2*_stats._kullbackLeibler[i][j])) / (_stats._nbSolidKmersPerBank[i]+_stats._nbSolidKmersPerBank[j]);

    return matrix;*/
//}

void SimkaDistance::get_abc(size_t bank1, size_t bank2, u_int64_t& a, u_int64_t& b, u_int64_t& c){

	a = _stats._matrixNbDistinctSharedKmers[bank1][bank2];
	b = (_stats._nbSolidDistinctKmersPerBank[bank1] - a);
	c = (_stats._nbSolidDistinctKmersPerBank[bank2] - a);

	//cout << bank1 << " " << bank2 << endl;
	//cout << "\t" << a << endl;
	//cout << "\t" << b << endl;
	//cout << "\t" << c << endl;
//void SimkaDistance::get_abc(SpeciesAbundanceVectorType& X, SpeciesAbundanceVectorType& Y, u_int64_t& a, u_int64_t& b, u_int64_t& c){
    /*
    gets classical values a, b and c as defined for diversity indices
    a: species preent in X and Y
    b: species present in X and not in Y
    c: species present in Y and not in X
    see Anderson  al., 2006, formulas (1) and (2)
    */
	/*
	cout << X.size() << endl;
	bitset<635956> pres_X;
	bitset<635956> pres_Y;
	//vector<bool> pres_X;
	//vector<bool> pres_Y;
	for(size_t i=0; i<X.size(); i++){
		pres_X.set(i, X[i] > 0);
		pres_Y.set(i, Y[i] > 0);
		//cout << (int)X[i] << endl;
	}
	//for(size_t i=0; i<X.size(); i++)

	a = (pres_X & pres_Y).count();
	pres_Y.flip();
	b = (pres_X & pres_Y).count();
	pres_Y.flip();
	pres_X.flip();
	c = (pres_Y & pres_X).count();

	cout << (a*100) / (double)X.size() << endl;
	cout << (b*100) / (double)X.size() << endl;
	cout << (c*100) / (double)X.size() << endl;
	*/
	/*
    pres_X = [abundance_2_presence(x) for x in X]
    pres_Y = [abundance_2_presence(y) for y in Y]
    x = np.array(pres_X)
    y = np.array(pres_Y)
    a = sum(1*np.logical_and(x,y))
    b = sum(1*np.logical_and(x,np.logical_not(y)))
    c = sum(1*np.logical_and(y,np.logical_not(x)))
    #
    return a, b, c*/
}

double SimkaDistance::brayCurtisSimilarity(size_t i, size_t j){
	/*
    """
    Bray Curtis index
    """
    n = len(X)
    bc_num = 0
    bc_den = 0
    for i in range(n):
        if (X[i] + Y[i] > 0):
            bc_num = bc_num + abs(X[i]-Y[i])
            bc_den = bc_den + X[i] + Y[i]
    bc = bc_num/bc_den
    #
    return bc*/
	/*
    u_int64_t bc_num = 0;
    u_int64_t bc_den = 0;
    for(size_t i=0; i<X.size(); i++){
        if(X[i] + Y[i] > 0){
            bc_num += abs(X[i]-Y[i]);
            bc_den += X[i] + Y[i];
        }
    }
    return bc_num/(double)bc_den;*/
	//cout <<  _stats._brayCurtisNumerator[i][j] << endl;
	//cout << _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j] << endl;
	double intersectionSize = _stats._brayCurtisNumerator[i][j];
	double unionSize = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j];

	/*
	cout << endl;
	cout << _stats._brayCurtisNumerator[i][j] << endl;
	cout << _stats._nbSolidKmersPerBank[i] << endl;
	cout << _stats._nbSolidKmersPerBank[j] << endl;
	cout << endl;*/

	return 100 * ((2*intersectionSize) / unionSize);
}

//double SimkaDistance::jaccardSimilarity(u_int64_t& a, u_int64_t& b, u_int64_t& c){
double SimkaDistance::jaccardSimilarity(size_t i, size_t j, SIMKA_MATRIX_TYPE type, SIMKA_PRESENCE_ABUNDANCE presenceOrAbundance){

	double intersectionSize = 0;
	double unionSize = 0;

	if(presenceOrAbundance == SIMKA_PRESENCE_ABUNDANCE::PRESENCE_ABSENCE){

	    if(type == SIMKA_MATRIX_TYPE::SYMETRICAL){
	    	intersectionSize = _stats._matrixNbDistinctSharedKmers[i][j];
	    	unionSize = _stats._nbSolidDistinctKmersPerBank[i] + _stats._nbSolidDistinctKmersPerBank[j] - intersectionSize;
	    }
	    else if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
	    	intersectionSize = _stats._matrixNbDistinctSharedKmers[i][j];
	    	unionSize = _stats._nbSolidDistinctKmersPerBank[i];
	    }
	}
	else if(presenceOrAbundance == SIMKA_PRESENCE_ABUNDANCE::ABUNDANCE){

	    if(type == SIMKA_MATRIX_TYPE::SYMETRICAL){
	    	intersectionSize = _stats._matrixNbSharedKmers[i][j] + _stats._matrixNbSharedKmers[j][i];
	    	unionSize = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j];
	    }
	    else if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
	    	intersectionSize = _stats._matrixNbSharedKmers[i][j];
	    	unionSize = _stats._nbSolidKmersPerBank[i];
	    }
	}


	return 100 * (intersectionSize / unionSize);

	/*
	double da = a;
	double db = b;
	double dc = c;

	return (db+dc) / (da+db+dc);

    """
    Jaccard index
    """
    a = float(abc[0])
    b = float(abc[1])
    c = float(abc[2])
    #
    jac = (b+c)/(a+b+c)
    #
    return jac*/
}

double SimkaDistance::sorensenSimilarity(u_int64_t& a, u_int64_t& b, u_int64_t& c, size_t i, size_t j, SIMKA_MATRIX_TYPE type){

	double intersectionSize = 0;
	double unionSize = 0;

    if(type == SIMKA_MATRIX_TYPE::SYMETRICAL){
    	intersectionSize = 2*a; //_stats._matrixNbDistinctSharedKmers[i][j];
    	unionSize = 2*a + b + c;//_stats._nbSolidDistinctKmersPerBank[i] + _stats._nbSolidDistinctKmersPerBank[j] - intersectionSize;
    }
    else if(type == SIMKA_MATRIX_TYPE::ASYMETRICAL){
    	intersectionSize = 2*a; // _stats._matrixNbDistinctSharedKmers[i][j];
    	unionSize = 2*a + b + c; //_stats._nbSolidDistinctKmersPerBank[i];
    }


	return 100 * ((intersectionSize) / unionSize);


	/*
	double da = a;
	double db = b;
	double dc = c;

	//return 2*da /
	return (db+dc) / (2*da+db+dc);

    """
    Jaccard index
    """
    a = float(abc[0])
    b = float(abc[1])
    c = float(abc[2])
    #
    jac = (b+c)/(a+b+c)
    #
    return jac*/
}

/*
double SimkaDistance::kullbackLeibler(SpeciesAbundanceVectorType& X, SpeciesAbundanceVectorType& Y){

	vector<double> XY(X.size());

    for(size_t i=0; i<X.size(); i++){
    	XY[i] = (X[i] + Y[i]) / 2.0;
    }

    double kl = 0;

    for(size_t i=0; i<X.size(); i++){
        if (X[i]*Y[i] > 0){
            double kl_x = X[i] * log(X[i]/XY[i]);
            double kl_y = Y[i] * log(Y[i]/XY[i]);
            kl += kl_x + kl_y;
        }
    }

    return kl;

}*/
