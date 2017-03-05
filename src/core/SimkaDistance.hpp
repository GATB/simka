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

#ifndef TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_
#define TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_

#include <gatb/gatb_core.hpp>
#include "../export/SimkaDistanceMatrixBinary.hpp"
#include "SimkaStatistics.hpp"

















class SimkaDistance {

public:

	vector<vector<float> > _matrix; //to remove
	vector<vector<float> > _matrix_rectangular;
	vector<vector<float> > _matrix_squaredHalf;
	SimkaStatistics& _stats;
	//SimkaDistanceParam _distanceParams;
	size_t _nbBanks;
	size_t _nbNewBanks;


	SimkaDistance(SimkaStatistics& stats) : _stats(stats){


		_nbBanks = _stats._nbBanks;
		_nbNewBanks = _stats._nbNewBanks;

		u_int64_t nbOldBanks = _nbBanks - _nbNewBanks;

		_matrix_rectangular.resize(_nbNewBanks);
		_matrix_squaredHalf.resize(_nbNewBanks-1);

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			_matrix_rectangular[i].resize(nbOldBanks);
		}

		u_int64_t size = _nbNewBanks-1;
		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			_matrix_squaredHalf[i].resize(size);
			size -= 1;
		}
	}

	void outputMatrix(const string& outputDirTemp, const vector<string>& bankNames, u_int64_t nbProcessedDatasets, u_int64_t nbNewDatasets){

		//SimkaDistance _simkaDistance(*this);


		SimkaDistanceMatrixBinary::saveMatrixIds(outputDirTemp, bankNames, nbProcessedDatasets, nbNewDatasets);


		compute_matrix_abundance_braycurtis();
		SimkaDistanceMatrixBinary::writeMatrixBinaryFromSplits(outputDirTemp, "mat_abundance_braycurtis", _matrix_rectangular, _matrix_squaredHalf);

		compute_matrix_presenceAbsence_jaccard();
		SimkaDistanceMatrixBinary::writeMatrixBinaryFromSplits(outputDirTemp, "mat_presenceAbsence_jaccard", _matrix_rectangular, _matrix_squaredHalf);

	}

	void clearMatrix(){
		/*
		for(size_t i=0; i<_matrix.size(); i++){
			for(size_t j=0; j<_matrix.size(); j++){
				_matrix[i][j] = 0;
			}
		}*/
	}


    void compute_matrix_abundance_braycurtis(){

    	clearMatrix();

		size_t nbOldBanks = _nbBanks - _nbNewBanks;

		//cout << "lala" << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular[0].size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf[0].size() << endl;

		for(size_t i=0; i<_stats._brayCurtisNumerator._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._brayCurtisNumerator._matrix_rectangular[i].size(); j++){

				double dist = distance_abundance_brayCurtis(i+nbOldBanks, j, i, j, _stats._brayCurtisNumerator._matrix_rectangular, _stats._nbKmersPerDatasetPairs._matrix_rectangular);
				_matrix_rectangular[i][j] = dist;
			}
		}

		for(size_t i=0; i<_stats._brayCurtisNumerator._matrix_squaredHalf.size(); i++){

			u_int64_t jOffset = _stats._brayCurtisNumerator._matrix_squaredHalf.size() - _stats._brayCurtisNumerator._matrix_squaredHalf[i].size();

			for(size_t j=0; j<_stats._brayCurtisNumerator._matrix_squaredHalf[i].size(); j++){

				double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._brayCurtisNumerator._matrix_squaredHalf, _stats._nbKmersPerDatasetPairs._matrix_squaredHalf);
				_matrix_squaredHalf[i][j] = dist;
			}
		}

    }

    void compute_matrix_presenceAbsence_jaccard(){
    	u_int64_t a, b, c;

    	clearMatrix();

		size_t nbOldBanks = _nbBanks - _nbNewBanks;

		//cout << "lala" << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular[0].size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf[0].size() << endl;

		for(size_t i=0; i<_stats._matrixNbDistinctSharedKmers._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._matrixNbDistinctSharedKmers._matrix_rectangular[i].size(); j++){


				//get_abc(i+nbOldBanks, j, i, j, _stats._matrixNbDistinctSharedKmers._matrix_rectangular, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(i+nbOldBanks, j, i, j, _stats._matrixNbDistinctSharedKmers._matrix_rectangular);

				//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j, i, j, _stats._brayCurtisNumerator._matrix_rectangular);
				_matrix_rectangular[i][j] = dist;
			}
		}

		for(size_t i=0; i<_stats._matrixNbDistinctSharedKmers._matrix_squaredHalf.size(); i++){

			u_int64_t jOffset = _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf.size() - _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf[i].size();

			for(size_t j=0; j<_stats._matrixNbDistinctSharedKmers._matrix_squaredHalf[i].size(); j++){

				//get_abc(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf, a, b ,c);
				double dist = distance_presenceAbsence_jaccardCanberra(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf);
				//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._brayCurtisNumerator._matrix_squaredHalf);
				_matrix_squaredHalf[i][j] = dist;
			}
		}

		/*
		//cout << _matrix.size() << " " << _matrix[0].size() << "   " << _nbNewBanks << endl;
		for(size_t i=0; i<nbOldBanks; i++){
			for(size_t j=0; j<_nbNewBanks; j++){

				//cout << i << " " << j << endl;
				get_abc(i, j, j+nbOldBanks, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);
				_matrix[i][j] = dist;

			}
		}

		for(size_t i=nbOldBanks; i<_nbBanks; i++){
			for(size_t j=(i-nbOldBanks)+1; j<_nbNewBanks; j++){

				get_abc(i, j, j+nbOldBanks, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);
				_matrix[i][j] = dist;


				size_t iTemp = i-nbOldBanks;
				size_t i2 = j;
				size_t j2 = iTemp;
				i2 += nbOldBanks;
				_matrix[i2][j2] = dist;

			}
		}*/



    }






private:


    //void get_abc(size_t i, size_t j, size_t i2, size_t j2, vector<vector<u_int64_t> >& crossedData, u_int64_t& a, u_int64_t& b, u_int64_t& c){

    	//a = crossedData[i2][j2];
    	//b = (_stats._nbSolidDistinctKmersPerBank[i] - a);
    	//c = (_stats._nbSolidDistinctKmersPerBank[j] - a);

    //}


    double distance_abundance_brayCurtis(size_t i, size_t j, size_t i2, size_t j2, vector<vector<u_int64_t> >& crossedData, vector<vector<u_int64_t> >& marginalData){
    	//if( i== 0 && j==1){
    	//	cout << crossedData[i2][j2] << endl;
    	//	cout << marginalData[i2][j2] << endl;
    	//}

    	double union_ = marginalData[i2][j2];
    	if(union_ == 0) return 1;

    	double intersection = 2 * crossedData[i2][j2];

    	double dist = 1 - intersection / union_;
    	/*
    	//double intersection = _stats._abundance_jaccard_intersection[i][j];
    	double union_ = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j];

    	if(union_ == 0) return 1;

    	double intersection = 2 * crossedData[i2][j2];

    	double jaccard = 1 - intersection / union_;
    	//double jaccard = intersection / union_;

    	//jaccard = (- 1 / (float)_stats._kmerSize) * log( ((float)(2*jaccard) / (float)(1+jaccard)));

    	//cout << jaccard << "   " << _stats._nbSolidKmersPerBank[i] << " " <<  _stats._nbSolidKmersPerBank[j] << "  " << _stats._brayCurtisNumerator[i][j] << endl;
    	return jaccard;*/
    	return dist;
    }


    double distance_presenceAbsence_jaccardCanberra(size_t i, size_t j, size_t i2, size_t j2, vector<vector<u_int64_t> >& crossedData){

    	//if( i== 0 && j==1){
    		//cout << crossedData[i2][j2] << endl;
    		//cout << _stats._matrixNbDistinctSharedKmers[i] << endl;
    		//cout << _stats._matrixNbDistinctSharedKmers[j] << endl;
    	//}
    	//double dist =
    	//a = crossedData[i2][j2];
    	//b = (_stats._nbSolidDistinctKmersPerBank[i] - a);
    	//c = (_stats._nbSolidDistinctKmersPerBank[j] - a);

    	/*
    	double a = (double) ua;
    	double b = (double) ub;
    	double c = (double) uc;

    	if((a+b+c) == 0) return 1;

    	double dist = (b+c) / (a+b+c);
    	//dist = 1 - dist;
    	//cout << _stats._kmerSize << endl;
    	//cout << (- 1 / (float)_stats._kmerSize) << "   " << ((float)(2*dist) / (float)(1+dist)) << "    " << log( ((float)(2*dist) / (float)(1+dist))) << endl;
    	//dist = (- 1 / (float)_stats._kmerSize) * log( ((float)(2*dist) / (float)(1+dist)));
    	return dist ;*/
    	return 1 - (crossedData[i2][j2] / (double) _stats._sketchSize);
    }






};




#endif /* TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_ */
