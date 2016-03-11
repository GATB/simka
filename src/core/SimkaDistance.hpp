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

const string STR_SIMKA_DISTANCE_BRAYCURTIS = "-bray-curtis";
const string STR_SIMKA_DISTANCE_CHORD = "-chord";
const string STR_SIMKA_DISTANCE_HELLINGER = "-hellinger";
const string STR_SIMKA_DISTANCE_CANBERRA = "-canberra";
const string STR_SIMKA_DISTANCE_KULCZYNSKI = "-kulczynski";


typedef vector<u_int16_t> SpeciesAbundanceVectorType;

enum SIMKA_MATRIX_TYPE{
	SYMETRICAL,
	ASYMETRICAL,
};

/*
class SimkaDistanceParam{

public:

	SimkaDistanceParam(){}

	SimkaDistanceParam(IProperties* params){
	    //_computeBrayCurtis = true;
		//_computeChord = true;
		//_computeHellinger = true;
		//_computeCanberra = true;
		//_computeKulczynski = true;
		//_computeBrayCurtis = params->get(STR_SIMKA_DISTANCE_BRAYCURTIS);
		//_computeChord = params->get(STR_SIMKA_DISTANCE_CHORD);
		//_computeHellinger = params->get(STR_SIMKA_DISTANCE_HELLINGER);
		//_computeCanberra = params->get(STR_SIMKA_DISTANCE_CANBERRA);
		//_computeKulczynski = params->get(STR_SIMKA_DISTANCE_KULCZYNSKI);
	}

	//bool _computeBrayCurtis;
	//bool _computeChord;
	//bool _computeHellinger;
	//bool _computeCanberra;
	//bool _computeKulczynski;
};*/


class SimkaStatistics{

public:

	SimkaStatistics(size_t nbBanks, bool computeSimpleDistances, bool computeComplexDistances);
	SimkaStatistics& operator+=  (const SimkaStatistics& other);
	void print();
	void load(const string& filename);
	void save(const string& filename);
	void outputMatrix(const string& outputDir, const vector<string>& _bankNames);

    size_t _nbBanks;
    bool _computeSimpleDistances;
    bool _computeComplexDistances;

	vector<u_int64_t> _nbSolidDistinctKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBank;

	vector<u_int64_t> _nbDistinctKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersSharedByBanksThreshold;

	vector<vector<u_int64_t> > _matrixNbDistinctSharedKmers;
	vector<vector<u_int64_t> > _matrixNbSharedKmers;

	vector<vector<u_int64_t> > _brayCurtisNumerator;
	//vector<vector<double> > _kullbackLeibler;


    //Abundance Chord
	vector<vector<long double> > _chord_NiNj;
	vector<long double> _chord_sqrt_N2;

    //Abundance Hellinger
	vector<vector<u_int64_t> > _hellinger_SqrtNiNj;
	vector<vector<u_int64_t> > _whittaker_minNiNj;
	vector<vector<long double> > _kullbackLeibler;

    //Abundance Canberra
	vector<vector<u_int64_t> > _canberra;

	//Abundance Kulczynski
	vector<vector<u_int64_t> > _kulczynski_minNiNj;


	//string _outputDir;

	u_int64_t _nbKmers;
	vector<u_int64_t> _nbKmersPerBank;
	u_int64_t _nbErroneousKmers;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSolidKmers;

	u_int64_t _nbSharedKmers;
	u_int64_t _nbDistinctSharedKmers;
	//SimkaDistanceParam _distanceParams;

	vector<u_int64_t> _datasetNbReads;
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

    vector<vector<float> > _matrixJaccardAbundance;
    vector<vector<float> > _matrixBrayCurtis;
    vector<vector<float> > _matrixChord;
    vector<vector<float> > _matrixHellinger;
    vector<vector<float> > _matrixWhittaker;
    vector<vector<float> > _matrixKullbackLeibler;
    vector<vector<float> > _matrixCanberra;
    vector<vector<float> > _matrixKulczynski;
    vector<vector<float> > _matrixSymJaccardAbundance;
    vector<vector<float> > _matrixAsymJaccardAbundance;
    vector<vector<float> > _matrixOchiai;
    vector<vector<float> > _matrixSorensen;

    //vector<vector<float> > _matrixKullbackLeibler;

    vector<vector<float> > _matrix_presenceAbsence_sorensenBrayCurtis;
    vector<vector<float> > _matrix_presenceAbsence_Whittaker;
    vector<vector<float> > _matrix_presenceAbsence_kulczynski;
    vector<vector<float> > _matrix_presenceAbsence_ochiai;
    vector<vector<float> > _matrix_presenceAbsence_chordHellinger;
    vector<vector<float> > _matrix_presenceAbsence_jaccardCanberra;
    vector<vector<float> > _matrix_presenceAbsence_jaccard_simka;
    vector<vector<float> > _matrix_presenceAbsence_jaccard_simka_asym;

private:


	vector<vector<float> > createSquaredMatrix(size_t n);
	void get_abc(size_t bank1, size_t bank2, u_int64_t& a, u_int64_t& b, u_int64_t& c);


    double distance_abundance_brayCurtis(size_t bank1, size_t bank2);
    double distance_abundance_chord(size_t i, size_t j);
    double distance_abundance_hellinger(size_t i, size_t j);
    double distance_abundance_whittaker(size_t i, size_t j);
    double distance_abundance_kullbackLeibler(size_t i, size_t j);
    double distance_abundance_canberra(size_t i, size_t j, u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_abundance_kulczynski(size_t i, size_t j);
    double distance_abundance_ochiai(size_t i, size_t j);
    double distance_abundance_sorensen(size_t i, size_t j);
    double distance_abundance_jaccard(size_t i, size_t j);
    double distance_abundance_jaccard_simka(size_t i, size_t j, SIMKA_MATRIX_TYPE type);


    double distance_presenceAbsence_chordHellinger(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_hellinger(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_whittaker(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_canberra(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_kulczynski(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_ochiai(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_sorensenBrayCurtis(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_jaccardCanberra(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_jaccard_simka(size_t i, size_t j, SIMKA_MATRIX_TYPE type);

	SimkaStatistics& _stats;
	//SimkaDistanceParam _distanceParams;
	size_t _nbBanks;

};

#endif /* TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_ */
