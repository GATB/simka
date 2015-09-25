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

#include "SimkaDistance.hpp"





SimkaStatistics::SimkaStatistics(size_t nbBanks, SimkaDistanceParam& distanceParams) :
_distanceParams(distanceParams)
{

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
	for(size_t i=0; i<_nbBanks; i++){
		_matrixNbDistinctSharedKmers[i].resize(nbBanks, 0);
		_matrixNbSharedKmers[i].resize(nbBanks, 0);
		_brayCurtisNumerator[i].resize(nbBanks, 0);
		//_kullbackLeibler[i].resize(nbBanks, 0);
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
    for(size_t i=0; i<_nbKmersSharedByBanksThreshold.size(); i++)
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

    for(size_t i=0; i<_nbBanks; i++){
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

void SimkaStatistics::load(Group& group){

    Storage::istream is (group, "simkaStats");

    //is.read ((char*)&_nbBanks,                sizeof(_nbBanks));
    is.read ((char*)&_nbKmers,                sizeof(_nbKmers));
    is.read ((char*)&_nbErroneousKmers,                sizeof(_nbErroneousKmers));
    is.read ((char*)&_nbDistinctKmers,                sizeof(_nbDistinctKmers));
    is.read ((char*)&_nbSolidKmers,                sizeof(_nbSolidKmers));

    is.read ((char*)_nbSolidDistinctKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    is.read ((char*)_nbKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    is.read ((char*)_nbSolidKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    is.read ((char*)_nbDistinctKmersSharedByBanksThreshold.data(), sizeof(u_int64_t)*_nbBanks);
    is.read ((char*)_nbKmersSharedByBanksThreshold.data(), sizeof(u_int64_t)*_nbBanks);

    for(size_t i=0; i<_nbBanks; i++){
    	is.read ((char*)_matrixNbDistinctSharedKmers[i].data(), sizeof(u_int64_t)*_nbBanks);
    	is.read ((char*)_matrixNbSharedKmers[i].data(), sizeof(u_int64_t)*_nbBanks);
    	is.read ((char*)_brayCurtisNumerator[i].data(), sizeof(u_int64_t)*_nbBanks);
    }

	/*
    tools::storage::impl::Storage::istream is (group, "simkaStats");

    is.read ((char*)&_nbpart,     sizeof(_nbpart));
    is.read ((char*)&_mm,         sizeof(_mm));
    is.read ((char*)&_nb_minims,  sizeof(_nb_minims));
    is.read ((char*)&_nbPass,     sizeof(_nbPass));

    DEBUG (("[Repartitor::load] :  _nbpart=%d  _mm=%d  _nb_minims=%d  _nbPass=%d \n",
        _nbpart, _mm, _nb_minims, _nbPass
    ));

    _repart_table.resize (_nb_minims);

    is.read ((char*)_repart_table.data(), sizeof(Value) * _nb_minims);

    is.read ((char*)&hasMinimizerFrequencies, sizeof(bool));

    u_int32_t magic = 0;
    is.read ((char*)&magic,  sizeof(magic));
    if (magic != MAGIC_NUMBER)  { throw system::Exception("Unable to load Repartitor (minimRepart), possibly due to bad format."); }

    if (hasMinimizerFrequencies)
    {
        tools::storage::impl::Storage::istream is2 (group, "minimFrequency");
        _freq_order = new uint32_t [_nb_minims];
        is2.read ((char*)_freq_order,     sizeof(uint32_t)*_nb_minims);

        is2.read ((char*)&magic,  sizeof(magic));
        if (magic != MAGIC_NUMBER)  { throw system::Exception("Unable to load Repartitor (minimFrequency), possibly due to bad format."); }
    }*/
}

void SimkaStatistics::save (Group& group){

    Storage::ostream os (group, "simkaStats");

    //os.write ((const char*)&_nbBanks,                sizeof(_nbBanks));
    os.write ((const char*)&_nbKmers,                sizeof(_nbKmers));
    os.write ((const char*)&_nbErroneousKmers,                sizeof(_nbErroneousKmers));
    os.write ((const char*)&_nbDistinctKmers,                sizeof(_nbDistinctKmers));
    os.write ((const char*)&_nbSolidKmers,                sizeof(_nbSolidKmers));

    os.write ((const char*)_nbSolidDistinctKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    os.write ((const char*)_nbKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    os.write ((const char*)_nbSolidKmersPerBank.data(), sizeof(u_int64_t)*_nbBanks);
    os.write ((const char*)_nbDistinctKmersSharedByBanksThreshold.data(), sizeof(u_int64_t)*_nbBanks);
    os.write ((const char*)_nbKmersSharedByBanksThreshold.data(), sizeof(u_int64_t)*_nbBanks);

    for(size_t i=0; i<_nbBanks; i++){
        os.write ((const char*)_matrixNbDistinctSharedKmers[i].data(), sizeof(u_int64_t)*_nbBanks);
        os.write ((const char*)_matrixNbSharedKmers[i].data(), sizeof(u_int64_t)*_nbBanks);
        os.write ((const char*)_brayCurtisNumerator[i].data(), sizeof(u_int64_t)*_nbBanks);
    }

    os.flush();
}

void SimkaStatistics::outputMatrix(const string& outputDir, const vector<string>& bankNames){


	SimkaDistance _simkaDistance(*this, _distanceParams);

	_outputFilenameSuffix = "";

	char buffer[200];

	//string strKmerSize = "_k";
	//snprintf(buffer,200,"%llu",_kmerSize);
	//strKmerSize += string(buffer);
	//_outputFilenameSuffix += strKmerSize;


	dumpMatrix(outputDir, bankNames, "mat_presenceAbsence_whittaker", _simkaDistance._matrix_presenceAbsence_Whittaker);
	dumpMatrix(outputDir, bankNames, "mat_presenceAbsence_kulczynski", _simkaDistance._matrix_presenceAbsence_kulczynski);
	dumpMatrix(outputDir, bankNames, "mat_presenceAbsence_sorensen", _simkaDistance._matrix_presenceAbsence_sorensen);

	dumpMatrix(outputDir, bankNames, "mat_presenceAbsence_jaccard", _simkaDistance._matrixSymJaccardPresenceAbsence);
	dumpMatrix(outputDir, bankNames, "mat_presenceAbsence_jaccard_asym", _simkaDistance._matrixAsymJaccardPresenceAbsence);

	dumpMatrix(outputDir, bankNames, "mat_abundance_jaccard", _simkaDistance._matrixSymJaccardAbundance);
	dumpMatrix(outputDir, bankNames, "mat_abundance_jaccard_asym", _simkaDistance._matrixAsymJaccardAbundance);

	if(_distanceParams._computeBrayCurtis)
		dumpMatrix(outputDir, bankNames, "mat_abundance_brayCurtis", _simkaDistance._matrixBrayCurtis);

	//dumpMatrix("mat_kullbackLeibler", _simkaDistance->getMatrixKullbackLeibler());

}



void SimkaStatistics::dumpMatrix(const string& outputDir, const vector<string>& bankNames, const string& outputFilename, const vector<vector<float> >& matrix){

	char buffer[200];
	string str;

	for(size_t i=0; i<matrix.size(); i++){
		str += ";" + bankNames[i];
		//str += ";" + datasetInfos[i]._name;
	}
	str += '\n';

	for(size_t i=0; i<matrix.size(); i++){

		str += bankNames[i] + ";";
		//str += datasetInfos[i]._name + ";";
		for(size_t j=0; j<matrix.size(); j++){

			//snprintf(buffer,200,"%.2f", matrix[i][j]);
			snprintf(buffer,200,"%f", matrix[i][j]);
			str += string(buffer) + ";";

			//str += to_string(matrix[i][j]) + ";";
		}

		//matrixNormalizedStr.erase(matrixNormalizedStr.end()-1);
		str.erase(str.size()-1);
		//str.pop_back(); //remove ; at the end of the line
		str += '\n';
	}


	gatb::core::system::IFile* file = gatb::core::system::impl::System::file().newFile(outputDir + "/" + outputFilename + _outputFilenameSuffix + ".csv", "wb");
	file->fwrite(str.c_str(), str.size(), 1);
	file->flush();
	delete file;

}



















SimkaDistance::SimkaDistance(SimkaStatistics& stats, SimkaDistanceParam& distanceParams) : _stats(stats), _distanceParams(distanceParams){
	_nbBanks = _stats._nbBanks;

	u_int64_t a;
	u_int64_t b;
	u_int64_t c;

    _matrixBrayCurtis = createSquaredMatrix(_nbBanks);
    _matrixSymJaccardPresenceAbsence = createSquaredMatrix(_nbBanks);
    _matrixAsymJaccardPresenceAbsence = createSquaredMatrix(_nbBanks);
    _matrixSymJaccardAbundance = createSquaredMatrix(_nbBanks);
    _matrixAsymJaccardAbundance = createSquaredMatrix(_nbBanks);

    _matrix_presenceAbsence_sorensen = createSquaredMatrix(_nbBanks);
    _matrix_presenceAbsence_Whittaker = createSquaredMatrix(_nbBanks);
    _matrix_presenceAbsence_kulczynski = createSquaredMatrix(_nbBanks);
    //_matrixAsymSorensen = createSquaredMatrix(_nbBanks);
    //_matrixKullbackLeibler = createSquaredMatrix(_nbBanks);

	for(size_t i=0; i<_nbBanks; i++){
		//SpeciesAbundanceVectorType& X_i = _stats._speciesAbundancePerDataset[i];

		//for(size_t j=0; j<_nbBanks; j++){
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

			//PresenceAbsence Sorensen
			double sorensen = distance_presenceAbsence_sorensen(a, b, c);
			_matrix_presenceAbsence_sorensen[i][j] = sorensen;
			_matrix_presenceAbsence_sorensen[j][i] = sorensen;


			//PresenceAbsence Whittaker
			double whittaker = distance_presenceAbsence_whittaker(a, b, c);
			_matrix_presenceAbsence_Whittaker[i][j] = whittaker;
			_matrix_presenceAbsence_Whittaker[j][i] = whittaker;


			//PresenceAbsence kulczynski
			double kulczynski = distance_presenceAbsence_kulczynski(a, b, c);
			_matrix_presenceAbsence_kulczynski[i][j] = kulczynski;
			_matrix_presenceAbsence_kulczynski[j][i] = kulczynski;




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

}




vector<vector<float> > SimkaDistance::createSquaredMatrix(size_t n){
    vector<vector<float> > matrix;

    matrix.resize(n);
    for(size_t i=0; i<n; i++)
    	matrix[i].resize(n, 100);

    return matrix;

}


void SimkaDistance::get_abc(size_t bank1, size_t bank2, u_int64_t& a, u_int64_t& b, u_int64_t& c){

	a = _stats._matrixNbDistinctSharedKmers[bank1][bank2];
	b = (_stats._nbSolidDistinctKmersPerBank[bank1] - a);
	c = (_stats._nbSolidDistinctKmersPerBank[bank2] - a);

}

double SimkaDistance::brayCurtisSimilarity(size_t i, size_t j){

	double intersectionSize = _stats._brayCurtisNumerator[i][j];
	double unionSize = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j];

	return 100 * ((2*intersectionSize) / unionSize);
}

//double SimkaDistance::jaccardSimilarity(u_int64_t& a, u_int64_t& b, u_int64_t& c){
double SimkaDistance::jaccardSimilarity(size_t i, size_t j, SIMKA_MATRIX_TYPE type, SIMKA_PRESENCE_ABUNDANCE presenceOrAbundance){

	double intersectionSize = 0;
	double unionSize = 0;

	if(presenceOrAbundance == PRESENCE_ABSENCE){

	    if(type == SYMETRICAL){
	    	intersectionSize = _stats._matrixNbDistinctSharedKmers[i][j];
	    	unionSize = _stats._nbSolidDistinctKmersPerBank[i] + _stats._nbSolidDistinctKmersPerBank[j] - intersectionSize;
	    }
	    else if(type == ASYMETRICAL){
	    	intersectionSize = _stats._matrixNbDistinctSharedKmers[i][j];
	    	unionSize = _stats._nbSolidDistinctKmersPerBank[i];
	    }
	}
	else if(presenceOrAbundance == ABUNDANCE){

	    if(type == SYMETRICAL){
	    	intersectionSize = _stats._matrixNbSharedKmers[i][j] + _stats._matrixNbSharedKmers[j][i];
	    	unionSize = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j];
	    }
	    else if(type == ASYMETRICAL){
	    	intersectionSize = _stats._matrixNbSharedKmers[i][j];
	    	unionSize = _stats._nbSolidKmersPerBank[i];
	    }
	}


	return 100 * (intersectionSize / unionSize);

}

double SimkaDistance::distance_presenceAbsence_sorensen(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc){

	double a = (double) ua;
	double b = (double) ub;
	double c = (double) uc;

	double distance = (b+c) / (2*a + b + c);

	return (1 - distance) * 100;
}

double SimkaDistance::distance_presenceAbsence_whittaker(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc){

	double a = (double) ua;
	double b = (double) ub;
	double c = (double) uc;

	double p1 = b / (a + b);
	double p2 = c / (a + c);

	double p3 = a / (a + b);
	double p4 = a / (a + c);

	double distance = 0.5 * (p1 + p2 + abs(p3-p4));

	return (1 - distance) * 100;
}

double SimkaDistance::distance_presenceAbsence_kulczynski(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc){

	double a = (double) ua;
	double b = (double) ub;
	double c = (double) uc;

	double p1 = a / (a + b);
	double p2 = a / (a + c);

	double distance = 1 - 0.5*(p1 + p2);

	return (1 - distance) * 100;
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
