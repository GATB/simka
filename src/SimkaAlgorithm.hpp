/*
 * SimkaAlgorithm.h
 *
 *  Created on: 18 mai 2015
 *      Author: gbenoit
 */

#ifndef TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_
#define TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_

#include <gatb/gatb_core.hpp>







/*********************************************************************
* ** SimkaCountProcessor
*********************************************************************/
template<size_t span>
class SimkaCountProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

	SimkaCountProcessor(size_t nbBanks, const pair<size_t, size_t>& abundanceThreshold);
	~SimkaCountProcessor();
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCountProcessor (_nbBanks, _abundanceThreshold);  }
	//CountProcessorAbstract<span>* clone ();
	void finishClones (vector<ICountProcessor<span>*>& clones);
	void finishClone(SimkaCountProcessor<span>* clone);
	virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum);

	bool isSolid(const CountVector& counts);
	void print();


	vector<u_int64_t> _nbSolidKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBankAbundance;

	vector<u_int64_t> _nbKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersAbundanceSharedByBanksThreshold;

	vector<vector<u_int64_t> > _matrixSharedKmers;
	vector<vector<u_int64_t> > _matrixSharedAbundanceKmers;

private:

    size_t         _nbBanks;
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
template<size_t span>
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
	void clear();

	string _outputDir;
	size_t _nbBanks;
	string _inputFilename;
	size_t _kmerSize;
	pair<size_t, size_t> _abundanceThreshold;
	//size_t _nbCores;

	string _banksInputFilename;
	IBank* _banks;
	IProperties* _options;

	SimkaCountProcessor<span>* _processor;
	vector<string> _bankNames;
};













#endif /* TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_ */
