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
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include<stdio.h>

//#define SIMKA_POTARA
//#define BOOTSTRAP
#define MAX_BOOTSTRAP 50
#define NB_BOOTSTRAP 45
//#define SIMKA_FUSION
//#define MULTI_PROCESSUS
//#define MULTI_DISK
//#define SIMKA_MIN
#include "SimkaDistance.hpp"

const string STR_SIMKA_SOLIDITY_PER_DATASET = "-solidity-single";
const string STR_SIMKA_MAX_READS = "-max-reads";
const string STR_SIMKA_MIN_READ_SIZE = "-min-read-size";
const string STR_SIMKA_MIN_READ_SHANNON_INDEX = "-read-shannon-index";
const string STR_SIMKA_MIN_KMER_SHANNON_INDEX = "-kmer-shannon-index";
const string STR_KMER_PER_READ = "-kmer-per-read";


enum SIMKA_SOLID_KIND{
	RANGE,
	SUM,
};



typedef u_int16_t bankIdType;

























class SimkaCounterBuilder
{
public:

    /** Constructor.
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
    SimkaCounterBuilder (size_t nbBanks=1)  :  _abundancePerBank(nbBanks)  {}

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank=0)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= 1;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank=0)  {  _abundancePerBank [idxBank] ++;  }

    /** Set the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void set (CountNumber val, size_t idxBank=0)  {  _abundancePerBank [idxBank] = val;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    CountNumber operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** */
    const CountVector& get () const { return _abundancePerBank; }

private:
    CountVector _abundancePerBank;
};




















/*********************************************************************
* ** SimkaCountProcessor
*********************************************************************/
template<size_t span>
class SimkaCountProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

	SimkaCountProcessor(SimkaStatistics& stats, size_t nbBanks, size_t kmerSize, const pair<size_t, size_t>& abundanceThreshold, SIMKA_SOLID_KIND solidKind, bool soliditySingle, double minKmerShannonIndex);
	~SimkaCountProcessor();
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCountProcessor (_stats, _nbBanks, _kmerSize, _abundanceThreshold, _solidKind, _soliditySingle, _minKmerShannonIndex);  }
	//CountProcessorAbstract<span>* clone ();
	void finishClones (vector<ICountProcessor<span>*>& clones);
	void finishClone(SimkaCountProcessor<span>* clone);
	virtual bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum);

	void computeStats(const CountVector& counts);
	//void updateBrayCurtis(int bank1, CountNumber abundance1, int bank2, CountNumber abundance2);

	inline bool isSolidVector(const CountVector& counts);
	//bool isSolid(CountNumber count);
	double getShannonIndex(const Type&  kmer);


private:

    size_t         _nbBanks;
    size_t _kmerSize;
	pair<CountNumber, CountNumber> _abundanceThreshold;
	bool isAbundanceThreshold;
    SIMKA_SOLID_KIND _solidKind;
    bool _soliditySingle;
    //IteratorListener* _progress;
    //vector<size_t> _countTotal;

	//u_int64_t _nbBanks;
    SimkaStatistics* _localStats;
    SimkaStatistics& _stats;
    u_int64_t _totalAbundance;

    u_int64_t _nbKmerCounted;
    double _minKmerShannonIndex;
    CountVector _solidCounts;

};

struct SimkaSequenceFilter
{
	u_int64_t _maxNbReads;
	//u_int64_t _maxNbReadsPerBank;
	u_int64_t _nbReadProcessed;
	//int* _bankIndex;
	//int* _datasetIndex;


	SimkaSequenceFilter(size_t minReadSize, double minShannonIndex){
		_maxNbReads = 0;
		_nbReadProcessed = 0;
		_minReadSize = minReadSize;
		_minShannonIndex = minShannonIndex;
	}

#ifdef BOOTSTRAP
	vector<bool> _bootstraps;


	void setBootstrap(vector<bool>& bootstraps){
		_bootstraps = bootstraps;
		//for(size_t i=0; i<_bootstraps.size(); i++)
		//	cout << _bootstraps[i];
		//cout << endl << endl;
	}

#endif

	void setMaxReads(u_int64_t maxReads){
		_maxNbReads = maxReads;
	}

	bool operator() (Sequence& seq){

		if(_maxNbReads != 0){
			if(_nbReadProcessed >= _maxNbReads){
				return false;
			}
		}

		//cout << seq.getIndex() << " " <<  _nbReadProcessed << endl;

#ifdef BOOTSTRAP
		int readPerBootstrap = _maxNbReads / MAX_BOOTSTRAP;
		int bootstrapIndex = seq.getIndex() / readPerBootstrap;
		if(!_bootstraps[bootstrapIndex]) return false;
		//cout << bootstrapIndex << endl;
#endif

		if(!isReadSizeValid(seq))
			return false;

		if(!isShannonIndexValid(seq))
			return false;


		//cout << _nbReadProcessed << endl;
		_nbReadProcessed += 1;
		return true;
	}

	bool isReadSizeValid(Sequence& seq){
		if(_minReadSize == 0) return true;
		return seq.getDataSize() >= _minReadSize;
	}

	bool isShannonIndexValid(Sequence& seq){
		if(_minShannonIndex == 0) return true;
		return getShannonIndex(seq) >= _minShannonIndex;
	}

	float getShannonIndex(Sequence& seq){

		static char nt2binTab[128] = {
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

	size_t _minReadSize;
	double _minShannonIndex;
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

	u_int64_t _numberRef;
	u_int64_t _totalSizeRef;
	u_int64_t _maxSizeRef;
    /** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
	SimkaBankFiltered (IBank* ref, const Filter& filter, const vector<u_int64_t>& nbReadsPerDataset, u_int64_t nbReadToProcess) : BankDelegate (ref), _filter(filter)  {

		_nbReadsPerDataset = nbReadsPerDataset;
		_nbReadToProcess = nbReadToProcess;

		ref->estimate(_numberRef, _totalSizeRef, _maxSizeRef);
	}

	/*
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize){

    	number = _nbReadToProcess;
    	totalSize = (_totalSizeRef*_nbReadToProcess)/_numberRef;
    	maxSize = _maxSizeRef;

    	//cout << number2 << endl;

    	//u_int64_t readSize = totalSize2 / number2;
    	//cout << "lal:" << number2 << endl;
    	//number = _maxReads;

    	//number = _nbReadToProcess;
    	//totalSize = _nbReadToProcess*readSize;
    	//maxSize = readSize;
    }*/

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

    	//cout << endl << "---" << endl;
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

            	//cout << _nbReadsPerDataset[i] << endl;
            	//Depending on the parameter -max-reads we truncate or not the reads iterator
            	if(_nbReadsPerDataset[i] == 0){

                	//Max nb reads parameter is set to 0. All the reads of each dataset are processed
                	iterators[i] = new FilterIterator<Sequence,Filter> (iterators[i], _filter);

            	}
            	else{

                	//We create a truncated iterator that stop processing reads when _nbReadsPerDataset[i] is reached
            		//cout << _nbReadsPerDataset[i] << endl;
            		SimkaTruncateIterator<Sequence>* truncIt = new SimkaTruncateIterator<Sequence>(iterators[i], _nbReadsPerDataset[i] + _nbReadsPerDataset[i]/5);
            		Filter filter(_filter);
            		filter.setMaxReads(_nbReadsPerDataset[i]);

#ifdef BOOTSTRAP

            		srand (time(NULL));
            		size_t nbBootstrap = 0;
            		vector<bool> iSBoostrap(MAX_BOOTSTRAP);

            		while(nbBootstrap != NB_BOOTSTRAP){
            			int index = rand() % iSBoostrap.size();

            			if(!iSBoostrap[index]){
            				iSBoostrap[index] = true;
            				nbBootstrap += 1;
            			}
            		}
            		filter.setBootstrap(iSBoostrap);

#endif
                	FilterIterator<Sequence,Filter>* filterIt = new FilterIterator<Sequence,Filter> (truncIt, filter);
                	iterators[i] = filterIt;

            	}

            }
            return new CompositeIterator<Sequence> (iterators);
        }
    }

private:

	vector<u_int64_t> _nbReadsPerDataset;
    Filter _filter;
    u_int64_t _nbReadToProcess;
};





/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class SimkaAlgorithm : public Algorithm
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    typedef typename ModelCanonical::Kmer                                   KmerType;

	SimkaAlgorithm(IProperties* options);
	~SimkaAlgorithm();
	void execute();
	void print();

	//void executeSimkamin();


    static string toString(u_int64_t value){
    	char buffer[40];
    	snprintf(buffer, 30, "%llu", value);
    	return string(buffer);
    }

private:

	void layoutInputFilename();
	void createBank();
	void count();

	void outputMatrix();

	//void dumpMatrix(const string& outputFilename, const vector<vector<float> >& matrix);
	//void outputHeatmap();
	//void __outputHeatmap(const string& outputFilenamePrefix, const string& matrixPercFilename, const string& matrixNormFilename);

	void clear();

	u_int64_t _maxMemory;
	size_t _nbCores;
	string _outputDir;
	string _outputDirTemp;
	size_t _nbBanks;
	string _inputFilename;
	size_t _kmerSize;
	pair<size_t, size_t> _abundanceThreshold;
	SIMKA_SOLID_KIND _solidKind;
	bool _soliditySingle;
	size_t _maxNbReads;
	size_t _minReadSize;
	double _minReadShannonIndex;
	double _minKmerShannonIndex;
	size_t _nbMinimizers;
	//size_t _nbCores;

	SimkaStatistics* _stats;
	//SimkaDistance* _simkaDistance;

	string _banksInputFilename;
	vector<string> _tempFilenamesToDelete;
	IBank* _banks;
	IProperties* _options;

	SimkaCountProcessor<span>* _processor;
	vector<string> _bankNames;
	vector<u_int64_t> _nbReadsPerDataset;

	string _outputFilenameSuffix;

	u_int64_t _totalKmers;
    vector<size_t> _nbBankPerDataset;
	//string _matDksNormFilename;
	//string _matDksPercFilename;
	//string _matAksNormFilename;
	//string _matAksPercFilename;
	//string _heatmapDksFilename;
	//string _heatmapAksFilename;

	/*
    gatb::core::tools::dp::IteratorListener* _progress;
    void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }


	size_t _nbPartitions;
    std::vector <std::vector<size_t> > _nbKmersPerPartitionPerBank;
    vector<vector<u_int64_t> > _nbk_per_radix_per_part;//number of kxmer per parti per rad


    Storage* _tmpPartitionsStorage;
    void setPartitionsStorage (Storage* tmpPartitionsStorage)  {  SP_SETATTR(tmpPartitionsStorage);  }

    Partition<Type>* _tmpPartitions;
    void setPartitions (Partition<Type>* tmpPartitions)  {  SP_SETATTR(tmpPartitions);  }

    vector<u_int64_t> _nbKmerPerPartitions;
    int getSizeofPerItem () const { return Type::getSize()/8 + sizeof(bankIdType); }
    std::vector<size_t> getNbCoresList();

    //this->_local_pInfo.incKmer_and_rad (p, radix_kxmer.getVal(), kx_size); //nb of superkmer per x per parti per radix

    //vector<SpeciesAbundanceVectorType > _speciesAbundancePerDataset;

    //MultiDiskStorage<Type>* _multiStorage;
    //u_int64_t _maxDisk;
	*/


};







#endif /* TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_ */
