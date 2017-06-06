/*
 * SimkaStatistics.hpp
 *
 *  Created on: 4 mars 2017
 *      Author: gbenoit
 */

#ifndef SIMKA2_SIMKA_SRC_CORE_SIMKASTATISTICS_HPP_
#define SIMKA2_SIMKA_SRC_CORE_SIMKASTATISTICS_HPP_

#include "DistanceMatrixData.hpp"


class SimkaStatistics{

public:

    size_t _nbBanks;
    size_t _nbNewBanks;
    size_t _symetricDistanceMatrixSize;
    size_t _kmerSize;
    size_t _sketchSize;

	DistanceMatrixData _brayCurtisNumerator;
	DistanceMatrixData _matrixNbDistinctSharedKmers;

	DistanceMatrixData _nbKmersPerDatasetPairs;
	vector<u_int64_t> _nbDistinctKmersPerDataset;
	//vector<u_int64_t> _nbKmersPerDataset;
	//DistanceMatrixData _nbKmersPerDatasetPairs;

	//vector<u_int64_t> _nbDistinctKmersPerDatasets;
	//vector<u_int64_t> _nbKmersPerDatasets;
    //vector<vector<u_int64_t> > _minHashUnionCross;
	//vector<vector<u_int64_t> > _brayCurtisNumerator;
	//vector<vector<double> > _kullbackLeibler;


	SimkaStatistics(size_t nbBanks, size_t nbNewBanks, bool computeSimpleDistances, bool computeComplexDistances, size_t kmerSize, size_t sketchSize)
	{

		_nbBanks = nbBanks;
		_nbNewBanks = nbNewBanks;
		_kmerSize = kmerSize;
		_sketchSize = sketchSize;

		_brayCurtisNumerator.resize(_nbBanks, _nbNewBanks);
		_matrixNbDistinctSharedKmers.resize(_nbBanks, _nbNewBanks);
		_nbKmersPerDatasetPairs.resize(_nbBanks, _nbNewBanks);
		_nbDistinctKmersPerDataset.resize(_nbBanks, 0);
		//_nbKmersPerDataset.resize(_nbBanks, 0);
		/*
		 *     	_minHashUnion.resize(_nbBanks, 0);
    	_minHashUnionCross.resize(_nbBanks);
    	for(size_t i=0; i<_minHashUnionCross.size(); i++){
    		_minHashUnionCross[i].resize(_nbBanks, 0);
    	}
		 */
	}


	SimkaStatistics& operator+=  (const SimkaStatistics& other){

		//cout << "lol1" << endl;

		size_t nbOldBanks = _nbBanks - _nbNewBanks;

		for(size_t i=0; i<_nbDistinctKmersPerDataset.size(); i++){
			_nbDistinctKmersPerDataset[i] += other._nbDistinctKmersPerDataset[i];
			//cout << i << " :" << _nbDistinctKmersPerDataset[i] << endl;
		}

		//cout << "lol2  " << _sketchSize << endl;
		//cout << "lala" << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular[0].size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf[0].size() << endl;

		//cout << _brayCurtisNumerator._matrix_rectangular.size() << endl;
		for(size_t i=0; i<_brayCurtisNumerator._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_brayCurtisNumerator._matrix_rectangular[i].size(); j++){

				size_t iCrossed = i;
				size_t jCrossed = j;
				size_t iMarginal = i+nbOldBanks;
				size_t jMarginal = j;

				//cout << "xd1  " << _nbDistinctKmersPerDataset[iMarginal] << " " << _nbDistinctKmersPerDataset[jMarginal] << " " << _matrixNbDistinctSharedKmers._matrix_squaredHalf[iCrossed][jCrossed] << endl;
				if(_nbDistinctKmersPerDataset[iMarginal] + _nbDistinctKmersPerDataset[jMarginal] - _matrixNbDistinctSharedKmers._matrix_rectangular[iCrossed][jCrossed] >= _sketchSize){
					continue;
				}
				//cout << "xd2" << endl;

				//cout << _brayCurtisNumerator._matrix_rectangular[i][j] << " " << _matrixNbDistinctSharedKmers._matrix_rectangular[i][j] << " " << _nbKmersPerDatasetPairs._matrix_rectangular[i][j] << endl;

				_brayCurtisNumerator._matrix_rectangular[i][j] += other._brayCurtisNumerator._matrix_rectangular[i][j];
				_matrixNbDistinctSharedKmers._matrix_rectangular[i][j] += other._matrixNbDistinctSharedKmers._matrix_rectangular[i][j];
				_nbKmersPerDatasetPairs._matrix_rectangular[i][j] += other._nbKmersPerDatasetPairs._matrix_rectangular[i][j];
				//cout << "xd3" << endl;
			}
		}

		//cout << "lol3" << endl;

		if(_brayCurtisNumerator._matrix_squaredHalf.size() > 0){
			for(size_t i=0; i<_brayCurtisNumerator._matrix_squaredHalf.size(); i++){

				u_int64_t jOffset = _brayCurtisNumerator._matrix_squaredHalf.size() - _brayCurtisNumerator._matrix_squaredHalf[i].size();

				for(size_t j=0; j<_brayCurtisNumerator._matrix_squaredHalf[i].size(); j++){


					size_t iCrossed = i;
					size_t jCrossed = j;
					size_t iMarginal = i+nbOldBanks;
					size_t jMarginal = j+nbOldBanks+1+jOffset;

					if(_nbDistinctKmersPerDataset[iMarginal] + _nbDistinctKmersPerDataset[jMarginal] - _matrixNbDistinctSharedKmers._matrix_squaredHalf[iCrossed][jCrossed] >= _sketchSize){
						continue;
					}


					_brayCurtisNumerator._matrix_squaredHalf[i][j] += other._brayCurtisNumerator._matrix_squaredHalf[i][j];
					_matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j] += other._matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j];
					_nbKmersPerDatasetPairs._matrix_squaredHalf[i][j] += other._nbKmersPerDatasetPairs._matrix_squaredHalf[i][j];

				}
			}
		}

		//cout << "lol4" << endl;
		/*
		for(size_t i=0; i<_brayCurtisNumerator._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_brayCurtisNumerator._matrix_rectangular[i].size(); j++){

				if(_nbDistinctKmersPerDataset[i] + _nbDistinctKmersPerDataset[j] - _matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j] >= _sketchSize){
					continue;
				}

				_brayCurtisNumerator._matrix_rectangular[i][j] += other._brayCurtisNumerator._matrix_rectangular[i][j];
				_matrixNbDistinctSharedKmers._matrix_rectangular[i][j] += other._matrixNbDistinctSharedKmers._matrix_rectangular[i][j];
				_nbKmersPerDatasetPairs._matrix_rectangular[i][j] += other._nbKmersPerDatasetPairs._matrix_rectangular[i][j];
			}
		}

		for(size_t i=0; i<_brayCurtisNumerator._matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_brayCurtisNumerator._matrix_squaredHalf[i].size(); j++){


				if(_nbDistinctKmersPerDataset[i] + _nbDistinctKmersPerDataset[j] - _matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j] >= _sketchSize){
					continue;
				}

				_brayCurtisNumerator._matrix_squaredHalf[i][j] += other._brayCurtisNumerator._matrix_squaredHalf[i][j];
				_matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j] += other._matrixNbDistinctSharedKmers._matrix_squaredHalf[i][j];
				_nbKmersPerDatasetPairs._matrix_squaredHalf[i][j] += other._nbKmersPerDatasetPairs._matrix_squaredHalf[i][j];
			}
		}
		*/

		//_nbKmersPerDatasetPairs += other._brayCurtisNumerator;
		//_matrixNbDistinctSharedKmers += other._matrixNbDistinctSharedKmers;
		//_nbKmersPerDatasetPairs += other._nbKmersPerDatasetPairs;

		return *this;
	}


	void load(const string& filename){

		IterableGzFile<long double>* file = new IterableGzFile<long double>(filename);
		Iterator<long double>* it = file->iterator();
		LOCAL(it);
		it->first();

		for(size_t i=0; i<_nbDistinctKmersPerDataset.size(); i++){
			_nbDistinctKmersPerDataset[i] = it->item();
			it->next();
		}

	    _brayCurtisNumerator.load(it);
	    _matrixNbDistinctSharedKmers.load(it);
	    _nbKmersPerDatasetPairs.load(it);



		delete file;

	}

	void save (const string& filename){


		BagGzFile<long double>* file = new BagGzFile<long double>(filename);

		for(size_t i=0; i<_nbDistinctKmersPerDataset.size(); i++){
			file->insert((long double)_nbDistinctKmersPerDataset[i]);
		}

	    _brayCurtisNumerator.save(file);
	    _matrixNbDistinctSharedKmers.save(file);
	    _nbKmersPerDatasetPairs.save(file);

		file->flush();

		delete file;

	}



};




#endif /* SIMKA2_SIMKA_SRC_CORE_SIMKASTATISTICS_HPP_ */
