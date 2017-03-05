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

		_brayCurtisNumerator += other._brayCurtisNumerator;
		_matrixNbDistinctSharedKmers += other._matrixNbDistinctSharedKmers;
		_nbKmersPerDatasetPairs += other._nbKmersPerDatasetPairs;

		return *this;
	}


	void load(const string& filename){

		IterableGzFile<long double>* file = new IterableGzFile<long double>(filename);
		Iterator<long double>* it = file->iterator();
		LOCAL(it);
		it->first();

	    _brayCurtisNumerator.load(it);
	    _matrixNbDistinctSharedKmers.load(it);
	    _nbKmersPerDatasetPairs.load(it);



		delete file;

	}

	void save (const string& filename){


		BagGzFile<long double>* file = new BagGzFile<long double>(filename);

	    _brayCurtisNumerator.save(file);
	    _matrixNbDistinctSharedKmers.save(file);
	    _nbKmersPerDatasetPairs.save(file);

		file->flush();

		delete file;

	}



};




#endif /* SIMKA2_SIMKA_SRC_CORE_SIMKASTATISTICS_HPP_ */
