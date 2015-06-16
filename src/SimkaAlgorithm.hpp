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

const string STR_SIMKA_SOLIDITY_PER_DATASET = "-solidity-single";
const string STR_SIMKA_MAX_READS = "-max-reads";
const string STR_SIMKA_MIN_READ_SIZE = "-min-read-size";
const string STR_SIMKA_MIN_SHANNON_INDEX = "-min-shannon-index";

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

struct SimkaSequenceFilter
{
	//u_int64_t _maxNbReads;
	//u_int64_t _maxNbReadsPerBank;
	//u_int64_t _nbReadProcessed;
	//int* _bankIndex;
	//int* _datasetIndex;

	SimkaSequenceFilter(size_t minReadSize, double minShannonIndex){

		_minReadSize = minReadSize;
		_minShannonIndex = minShannonIndex;
	}

	bool operator() (Sequence& seq){

		if(!isReadSizeValid(seq))
			return false;

		if(!isShannonIndexValid(seq))
			return false;

		return true;
	}

	bool isReadSizeValid(Sequence& seq){
		if(_minReadSize == 0) return true;
		return seq.getDataSize() >= _minReadSize;
	}

	float getShannonIndex(Sequence& seq){
		float index = 0;
		//float freq [5];

		vector<float> _freqs(5, 0);

		char* seqStr = seq.getDataBuffer();

		// Frequency of each letter (A, C, G, T or N)
		for(size_t i=0; i < seq.getDataSize(); i++)
			_freqs[nt2binTab[(unsigned char)seqStr[i]]] += 1.0;

		// Shannon index calculation
		for (size_t i=0; i<_freqs.size(); i++){
			_freqs[i] /= (float) seq.getDataSize();
			if (_freqs[i] != 0)
				index += _freqs[i] * log (_freqs[i]) / log(2);
		}
		return abs(index);

	}

	bool isShannonIndexValid(Sequence& seq){
		if(_minShannonIndex == 0) return true;
		return getShannonIndex(seq) >= _minShannonIndex;
	}

	size_t _minReadSize;
	double _minShannonIndex;
	const char nt2binTab[128] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0, 0, //69
		0, 3, 0, 0, 0, 0, 0, 0, 4, 0, //79
		0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		};
};


template <class Item> class SimkaTruncateIterator : public TruncateIterator<Item>
{
public:

	SimkaTruncateIterator (Iterator<Item>* ref, u_int64_t limit, bool initRef=true)
        : TruncateIterator<Item>(*ref, limit, initRef), _ref2(0){ setRef(ref); }

private:

    Iterator<Item>* _ref2;
    void setRef (Iterator<Item>* ref2)  { SP_SETATTR(ref2); }

};

template<typename Filter> class SimkaBankFiltered : public BankDelegate
{
public:

    /** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
	SimkaBankFiltered (IBank* ref, const Filter& filter, const vector<u_int64_t>& nbReadsPerDataset) : BankDelegate (ref), _filter(filter)  {
		_nbReadsPerDataset = nbReadsPerDataset;
	}

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

    	//cout << "lala" << endl;
        // We create one iterator from the reference
        Iterator<Sequence>* it = _ref->iterator ();

        // We get the composition for this iterator
        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

        if (iterators.size() == 1)  { return new FilterIterator<Sequence,Filter> (it, _filter); }
        else
        {
            // We are going to create a new CompositeIterator, we won't need the one we just got from the reference
            LOCAL(it);

            // We may have to encapsulate each sub iterator with the filter.
            for (size_t i=0; i<iterators.size(); i++)  {

            	//cout << "\t\t" << _nbReadsPerDataset[i] << endl;

            	//Depending on the parameter -max-reads we truncate or not the reads iterator
            	if(_nbReadsPerDataset[i] == 0){

                	//Max nb reads parameter is set to 0. All the reads of each dataset are processed
                	iterators[i] = new FilterIterator<Sequence,Filter> (iterators[i], _filter); ;

            	}
            	else{

                	//We create a truncated iterator that stop processing reads when _nbReadsPerDataset[i] is reached
            		SimkaTruncateIterator<Sequence>* truncIt = new SimkaTruncateIterator<Sequence>(iterators[i], _nbReadsPerDataset[i]);
                	FilterIterator<Sequence,Filter>* filterIt = new FilterIterator<Sequence,Filter> (truncIt, _filter);
                	iterators[i] = filterIt;

            	}

            }
            return new CompositeIterator<Sequence> (iterators);
        }
    }

private:

	vector<u_int64_t> _nbReadsPerDataset;
    Filter _filter;
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
	void createBank();
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
	size_t _maxNbReads;
	size_t _minReadSize;
	double _minShannonIndex;
	//size_t _nbCores;

	int _datasetIndex;
	int _bankIndex;

	string _banksInputFilename;
	vector<string> _tempFilenamesToDelete;
	IBank* _banks;
	IProperties* _options;

	SimkaCountProcessor<span>* _processor;
	vector<string> _bankNames;
	vector<u_int64_t> _nbReadsPerDataset;

	string _matDksNormFilename;
	string _matDksPercFilename;
	string _matAksNormFilename;
	string _matAksPercFilename;
	string _heatmapDksFilename;
	string _heatmapAksFilename;



};













#endif /* TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_ */
