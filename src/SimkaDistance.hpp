/*
 * SimkaDistance.hpp
 *
 *  Created on: 17 juin 2015
 *      Author: gbenoit
 */

#ifndef TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_
#define TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_

#include <gatb/gatb_core.hpp>

enum SIMKA_MATRIX_TYPE{
	NORMALIZED,
	ASYMETRICAL
};

class SimkaStatistics{

public:

	SimkaStatistics(size_t nbBanks);
	SimkaStatistics& operator+=  (const SimkaStatistics& other);
	void print();

    size_t _nbBanks;

	vector<u_int64_t> _nbSolidDistinctKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBank;

	vector<u_int64_t> _nbDistinctKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersSharedByBanksThreshold;

	vector<vector<u_int64_t> > _matrixNbDistinctSharedKmers;
	vector<vector<u_int64_t> > _matrixNbSharedKmers;

	vector<vector<u_int64_t> > _brayCurtisNumerator;
	vector<vector<double> > _kullbackLeibler;


	//string _outputDir;

	u_int64_t _nbKmers;
	vector<u_int64_t> _nbKmersPerBank;
	u_int64_t _nbErroneousKmers;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSolidKmers;

	//u_int64_t _nbKmersInCoupleBankSupRatio;

	//unordered_map<string, histo_t> _histos;
};


class SimkaDistance {

public:

	SimkaDistance(SimkaStatistics& stats);
	//virtual ~SimkaDistance();

	vector<vector<float> > getMatrixSorensen(SIMKA_MATRIX_TYPE type);
	vector<vector<float> > getMatrixJaccard();
	vector<vector<float> > getMatrixAKS(SIMKA_MATRIX_TYPE type);
	vector<vector<float> > getMatrixBrayCurtis();
	vector<vector<float> > getMatrixKullbackLeibler();


private:


	vector<vector<float> > createSquaredMatrix(int n);
	//void outputMatrix();


	SimkaStatistics& _stats;
	size_t _nbBanks;

};

#endif /* TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_ */
