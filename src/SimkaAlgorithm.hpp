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

#ifndef TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_
#define TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_

#include <gatb/gatb_core.hpp>
#include<stdio.h>

const string STR_SOLIDITY_PER_DATASET = "-solidity-single";

enum SIMKA_SOLID_KIND{
	RANGE,
	SUM,
};



/*********************************************************************
* ** SimkaCountProcessor
*********************************************************************/
template<size_t span>
class SimkaCountProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

	SimkaCountProcessor(size_t nbBanks, const pair<size_t, size_t>& abundanceThreshold, SIMKA_SOLID_KIND solidKind, bool soliditySingle);
	~SimkaCountProcessor();
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCountProcessor (_nbBanks, _abundanceThreshold, _solidKind, _soliditySingle);  }
	//CountProcessorAbstract<span>* clone ();
	void finishClones (vector<ICountProcessor<span>*>& clones);
	void finishClone(SimkaCountProcessor<span>* clone);
	virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum);
	void computeStats(const CountVector& counts);

	bool isSolidVector(const CountVector& counts);
	bool isSolid(CountNumber count);
	void print();

	vector<u_int64_t> _nbSolidKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBankAbundance;

	vector<u_int64_t> _nbKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersAbundanceSharedByBanksThreshold;

	vector<vector<u_int64_t> > _matrixSharedKmers;
	vector<vector<u_int64_t> > _matrixSharedAbundanceKmers;

private:

    size_t         _nbBanks;
    SIMKA_SOLID_KIND _solidKind;
    bool _soliditySingle;
    //vector<size_t> _countTotal;

	//u_int64_t _nbBanks;
	pair<size_t, size_t> _abundanceThreshold;
	//string _outputDir;

	u_int64_t _nbKmers;
	vector<u_int64_t> _nbKmersPerBank;
	u_int64_t _nbErroneousKmers;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSolidKmers;

	//u_int64_t _nbKmersInCoupleBankSupRatio;

	//unordered_map<string, histo_t> _histos;

};




/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class SimkaAlgorithm {

public:

	SimkaAlgorithm(IProperties* options);
	~SimkaAlgorithm();
	void execute();


private:

	void layoutInputFilename();
	void count();
	void outputMatrix();
	void dumpMatrix(const string& outputFilename, vector<vector<float> >& matrix);
	void outputHeatmap();
	void __outputHeatmap(const string& matrixPercFilename, const string& matrixNormFilename);
	void printHelp();
	void clear();

	string _outputDir;
	size_t _nbBanks;
	string _inputFilename;
	size_t _kmerSize;
	pair<size_t, size_t> _abundanceThreshold;
	SIMKA_SOLID_KIND _solidKind;
	bool _soliditySingle;
	//size_t _nbCores;

	string _banksInputFilename;
	vector<string> _tempFilenamesToDelete;
	IBank* _banks;
	IProperties* _options;

	SimkaCountProcessor<span>* _processor;
	vector<string> _bankNames;

	string _matDksNormFilename;
	string _matDksPercFilename;
	string _matAksNormFilename;
	string _matAksPercFilename;
	string _heatmapDksFilename;
	string _heatmapAksFilename;

};













#endif /* TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_ */
