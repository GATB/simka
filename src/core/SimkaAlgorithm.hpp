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
#include "SimkaUtils.hpp"

//#define PRINT_STATS
//#define CHI2_TEST
//#define SIMKA_POTARA
//#define BOOTSTRAP
#define MAX_BOOTSTRAP 50
#define NB_BOOTSTRAP 45
//#define SIMKA_FUSION
//#define MULTI_PROCESSUS
//#define MULTI_DISK
//#define SIMKA_MIN
#include "SimkaDistance.hpp"



























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

























/*
template <class Item> class SimkaTruncateIterator : public TruncateIterator<Item>
{
public:

	SimkaTruncateIterator (Iterator<Item>* ref, u_int64_t limit, bool initRef=true)
        : TruncateIterator<Item>(*ref, limit, initRef), _ref2(0){ setRef(ref); }

private:

    Iterator<Item>* _ref2;
    void setRef (Iterator<Item>* ref2)  { SP_SETATTR(ref2); }

};*/

template<typename Filter> class SimkaBankFiltered : public BankDelegate
{
public:

	u_int64_t _refNbReads;
	u_int64_t _refTotalSeqSize;
	u_int64_t _refMaxReadSize;

	/** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
	SimkaBankFiltered (IBank* ref, const Filter& filter, const vector<size_t>& nbPaireds, u_int64_t maxReads) : BankDelegate (ref), _filter(filter)  {

		_nbPaireds = nbPaireds;
		_maxReads = maxReads;
		_nbBanks = ref->getCompositionNb();
		ref->estimate(_refNbReads, _refTotalSeqSize, _refMaxReadSize);


    	//cout << _refNbReads << endl;
		//cout << _refTotalSeqSize << endl;
		//cout << _refMaxReadSize << endl;
	}


    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize){


    	if(_maxReads == 0){
    		number = _refNbReads;
    		totalSize = _refTotalSeqSize;
    		maxSize = _refMaxReadSize;
    	}
    	else{

    		u_int64_t maxReads = 0;
    		for(size_t i=0; i<_nbBanks; i++){
    			maxReads += _maxReads * _nbPaireds[i];
    		}
    		//cout << _refNbReads << endl;
    		//cout << _maxReads*_nbBanks << endl;
    		maxReads = min (maxReads, _refNbReads);
			//cout << "ha " <<  maxReads << endl;

			if(maxReads == _refNbReads){
	    		number = _refNbReads;
	    		totalSize = _refTotalSeqSize;
	    		maxSize = _refMaxReadSize;
			}
			else{
				number = maxReads;
				double factor =  (double)maxReads / (double)_refNbReads;
				totalSize = _refTotalSeqSize * factor;
				maxSize = _refMaxReadSize;
			}
    	}

    	//number = _maxReads;
    	//totalSize = (_totalSizeRef*_nbReadToProcess)/_numberRef;
    	//maxSize = _maxSizeRef;

    	//cout << number2 << endl;

    	//u_int64_t readSize = totalSize2 / number2;
    	//cout << "lal:" << number2 << endl;
    	//number = _maxReads;

    	//number = _nbReadToProcess;
    	//totalSize = _nbReadToProcess*readSize;
    	//maxSize = readSize;

    	cout << number << endl;
    	//cout << totalSize << endl;
    	//cout << maxSize << endl;
    }

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

    	//cout << endl << "---" << endl;
    	//cout << "lala" << endl;
        // We create one iterator from the reference
        Iterator<Sequence>* it = _ref->iterator ();

        // We get the composition for this iterator
        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

        //if (iterators.size() == 1)  { return new FilterIterator<Sequence,Filter> (it, _filter); }
        //else
        //{
            // We are going to create a new CompositeIterator, we won't need the one we just got from the reference
		LOCAL(it);

		// We may have to encapsulate each sub iterator with the filter.
		for (size_t i=0; i<iterators.size(); i++)  {

			/*
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

				//CancellableIterator<Sequence>* truncIt = new CancellableIterator<Sequence>(*iterators[i]);
				Filter filter(_filter);
				//filter.setMaxReads(_nbReadsPerDataset[i]);
				//filter.setIt(truncIt);

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
				FilterIterator<Sequence,Filter>* filterIt = new FilterIterator<Sequence,Filter> (iterators[i], filter);
				iterators[i] = filterIt;

			}*/


				//Iterator<Sequence>* it = iterators[i];
				//std::vector<Iterator<Sequence>*> iterators_ = it->getComposition();
				iterators[i] = new SimkaInputIterator<Sequence, Filter> (iterators[i], _nbPaireds[i], _maxReads, _filter);

            }

		return new CompositeIterator<Sequence> (iterators);
    }

private:

	vector<size_t> _nbPaireds;
    Filter _filter;
    u_int64_t _maxReads;
    size_t _nbBanks;
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

protected:


    bool setup();
    bool isInputValid();
    void parseArgs();
    bool createDirs();
    void computeMaxReads();
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
	pair<CountNumber, CountNumber> _abundanceThreshold;
	SIMKA_SOLID_KIND _solidKind;
	bool _soliditySingle;
	int64_t _maxNbReads;
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

	vector<string> _bankNames;
	//vector<u_int64_t> _nbReadsPerDataset;

	string _outputFilenameSuffix;

	u_int64_t _totalKmers;
    vector<size_t> _nbBankPerDataset;

	string _largerBankId;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	bool _keepTmpFiles;
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
