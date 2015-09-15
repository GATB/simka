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

typedef vector<u_int16_t> SpeciesAbundanceVectorType;

enum SIMKA_MATRIX_TYPE{
	SYMETRICAL,
	ASYMETRICAL,
};

enum SIMKA_PRESENCE_ABUNDANCE{
	PRESENCE_ABSENCE,
	ABUNDANCE,
};

class SimkaStatistics{

public:

	SimkaStatistics(size_t nbBanks);
	SimkaStatistics& operator+=  (const SimkaStatistics& other);
	void print();
	void load(Group& group);
	void save(Group& group);
	void outputMatrix(const string& outputDir, const vector<string>& _bankNames);

    size_t _nbBanks;

	vector<u_int64_t> _nbSolidDistinctKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBank;

	vector<u_int64_t> _nbDistinctKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersSharedByBanksThreshold;

	vector<vector<u_int64_t> > _matrixNbDistinctSharedKmers;
	vector<vector<u_int64_t> > _matrixNbSharedKmers;

	vector<vector<u_int64_t> > _brayCurtisNumerator;
	//vector<vector<double> > _kullbackLeibler;


	//string _outputDir;

	u_int64_t _nbKmers;
	vector<u_int64_t> _nbKmersPerBank;
	u_int64_t _nbErroneousKmers;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSolidKmers;

	//u_int64_t _nbKmersInCoupleBankSupRatio;

	//unordered_map<string, histo_t> _histos;
private:

	void dumpMatrix(const string& outputDir, const vector<string>& _bankNames, const string& outputFilename, const vector<vector<float> >& matrix);
	string _outputFilenameSuffix;
};


class SimkaDistance {

public:

	SimkaDistance(SimkaStatistics& stats);
	//virtual ~SimkaDistance();

	//vector<vector<float> > getMatrixSorensen(SIMKA_MATRIX_TYPE type);
	//vector<vector<float> > getMatrixJaccard();
	//vector<vector<float> > getMatrixAKS(SIMKA_MATRIX_TYPE type);
	//vector<vector<float> > getMatrixBrayCurtis();
	//vector<vector<float> > getMatrixKullbackLeibler();

    vector<vector<float> > _matrixSymJaccardPresenceAbsence;
    vector<vector<float> > _matrixAsymJaccardPresenceAbsence;
    vector<vector<float> > _matrixSymJaccardAbundance;
    vector<vector<float> > _matrixAsymJaccardAbundance;
    vector<vector<float> > _matrixSymSorensen;
    vector<vector<float> > _matrixAsymSorensen;
    vector<vector<float> > _matrixBrayCurtis;
    //vector<vector<float> > _matrixKullbackLeibler;

private:


	vector<vector<float> > createSquaredMatrix(size_t n);
	void get_abc(size_t bank1, size_t bank2, u_int64_t& a, u_int64_t& b, u_int64_t& c);
    //void get_abc(SpeciesAbundanceVectorType& X, SpeciesAbundanceVectorType& Y, u_int64_t& a, u_int64_t& b, u_int64_t& c);
    double jaccardSimilarity(size_t i, size_t j, SIMKA_MATRIX_TYPE type, SIMKA_PRESENCE_ABUNDANCE presenceOrAbundance);
    double sorensenSimilarity(u_int64_t& a, u_int64_t& b, u_int64_t& c, size_t i, size_t j, SIMKA_MATRIX_TYPE type);
    double brayCurtisSimilarity(size_t bank1, size_t bank2);
    //double kullbackLeibler(SpeciesAbundanceVectorType& X, SpeciesAbundanceVectorType& Y);
	//void outputMatrix();


	SimkaStatistics& _stats;
	size_t _nbBanks;

};

#endif /* TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_ */
