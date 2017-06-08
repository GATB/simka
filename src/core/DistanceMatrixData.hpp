/*
 * SimkaDistanceMatrix.h
 *
 *  Created on: 4 mars 2017
 *      Author: gbenoit
 */

#ifndef SIMKA2_SIMKA_SRC_CORE_DISTANCEMATRIXDATA_HPP_
#define SIMKA2_SIMKA_SRC_CORE_DISTANCEMATRIXDATA_HPP_

//#include <vector>


class DistanceMatrixData
{

public:

	//size_t _nbBanks;
	//size_t _nbNewBanks;

	vector<vector<u_int64_t> > _matrix_rectangular;
	vector<vector<u_int64_t> > _matrix_squaredHalf;

	DistanceMatrixData(){};

	void resize(u_int64_t nbBanks, u_int64_t nbNewBanks){

		//cout << "Create distance matrix data" << endl;
		//cout << "\t" << nbBanks << " " << nbNewBanks << endl << endl;
		//_nbBanks = nbBanks;
		//_nbNewBanks = nbNewBanks;

		u_int64_t nbOldBanks = nbBanks - nbNewBanks;

		if(nbOldBanks > 0){
			_matrix_rectangular.resize(nbNewBanks);
		}
		else{
			_matrix_rectangular.resize(0);
		}

		_matrix_squaredHalf.resize(nbNewBanks-1);

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			_matrix_rectangular[i].resize(nbOldBanks);
		}

		u_int64_t size = nbNewBanks-1;
		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			_matrix_squaredHalf[i].resize(size);
			size -= 1;
		}
	}

	DistanceMatrixData& operator+=  (const DistanceMatrixData& other){
		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				_matrix_rectangular[i][j] += other._matrix_rectangular[i][j];
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				_matrix_squaredHalf[i][j] += other._matrix_squaredHalf[i][j];
			}
		}

		return *this;
	}

	void load(Iterator<u_int64_t>* gzIt){
		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				_matrix_rectangular[i][j] = gzIt->item();
				gzIt->next();
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				_matrix_squaredHalf[i][j] = gzIt->item();
				gzIt->next();
			}
		}
	}

	void save (Bag<u_int64_t>* bag){

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				bag->insert((long double)_matrix_rectangular[i][j]);
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				bag->insert((long double)_matrix_squaredHalf[i][j]);
			}
		}

	}

};


#endif /* SIMKA2_SIMKA_SRC_CORE_DISTANCEMATRIXDATA_HPP_ */
