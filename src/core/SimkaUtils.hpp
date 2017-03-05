/*
 * Utils.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_CORE_SIMKAUTILS_HPP_
#define GATB_SIMKA_SRC_CORE_SIMKAUTILS_HPP_

#include <gatb/gatb_core.hpp>
#include "SimkaDistance.hpp"

//#define _sketchSize 1000000

const string STR_SIMKA_SOLIDITY_PER_DATASET = "-solidity-single";
const string STR_SIMKA_MAX_READS = "-max-reads";
const string STR_SIMKA_MIN_READ_SIZE = "-min-read-size";
const string STR_SIMKA_MIN_READ_SHANNON_INDEX = "-read-shannon-index";
const string STR_SIMKA_MIN_KMER_SHANNON_INDEX = "-kmer-shannon-index";
const string STR_KMER_PER_READ = "-kmer-per-read";
const string STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES= "-simple-dist";
const string STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES = "-complex-dist";
const string STR_SIMKA_KEEP_TMP_FILES = "-keep-tmp";
const string STR_SIMKA_COMPUTE_DATA_INFO = "-data-info";

#define SIMKA2_LZ4_CACHE_NB_ITEMS 1000

enum SIMKA_SOLID_KIND{
	RANGE,
	SUM,
};



typedef u_int16_t bankIdType;



class SimkaUtils {
public:
	SimkaUtils(){

	}

	~SimkaUtils(){

	}
};



/********************************************************************************/
/**
 *
 */
template <class Item, typename Filter> class SimkaInputIterator : public Iterator<Item>
{
public:

	/** Constructor.
	* \param[in] ref : the referred iterator
	* \param[in] initRef : will call 'first' on the reference if true
	*/
	SimkaInputIterator(Iterator<Item>* refs, size_t nbBanks, u_int64_t maxReads, Filter filter)
	:  _filter(filter), _mainref(0) {

		setMainref(refs);
		_ref = _mainref->getComposition()[0];
		_isDone = false;
		_nbDatasets = nbBanks;
		_nbBanks = _mainref->getComposition().size() / _nbDatasets;
		_maxReads = maxReads;
		_nbReadProcessed = 0;
		_currentBank = 0;
		_currentInternalBank = 0;
		_currentDataset = 0;

	}


    bool isFinished(){
        if(_currentDataset == _nbDatasets){
                _isDone = true;
                return true;
        }
        return false;
    }

	void nextDataset(){
		_currentDataset += 1;

		if(isFinished()) return;

		_currentBank = _currentDataset * _nbBanks;

		_currentInternalBank = 0;
		_nbReadProcessed = 0;

		if(isFinished()) return;

		_ref = _mainref->getComposition()[_currentBank];
		_isDone = false;
		first();
		//nextBank();
	}

	void nextBank(){
		//cout << "next bank" << endl;
		//cout << "next bank "<< endl;
		_currentInternalBank += 1;
		if(_currentInternalBank == _nbBanks){
			nextDataset();
		}
		else{
			_isDone = false;
			_currentBank += 1;
			_ref = _mainref->getComposition()[_currentBank];
			first();
		}
	}

    void first()
    {

        _ref->first();

        while (!_ref->isDone() && _filter(_ref->item())==false)
                _ref->next();

        _isDone = _ref->isDone();

        if(!_isDone) *(this->_item) = _ref->item();

    }

	void next(){


		if(isFinished()){
			_isDone = true;
			return;
		}

		//cout << "haha" << endl;

		_ref->next();
		while (!_ref->isDone() && _filter(_ref->item())==false) _ref->next();

		_isDone = _ref->isDone();

		//cout << "haha" << endl;
		//if(!_isDone){
			//cout << _currentBank << "  " << _isDone << endl;

		//}

		//cout << _nbReadProcessed << "  " << _currentBank << "    " << _nbBanks << "   " << _maxReads << endl;


		if(_isDone){
			if(isFinished()){
				//cout << _nbReadProcessed << endl;
				return;
			}
			else{
				//cout << _nbReadProcessed << endl;
				nextBank();
				if(isFinished()){
					//cout << _nbReadProcessed << endl;
					return;
				}
			}
		}
		else{
			*(this->_item) = _ref->item();
			_nbReadProcessed += 1;
		}

		if(_maxReads && _nbReadProcessed >= _maxReads){
			if(isFinished())
				return;
			else
				nextDataset();
		}

	}

    /** \copydoc  Iterator::isDone */
    bool isDone()  {  return _isDone;  }

    /** \copydoc  Iterator::item */
    Item& item ()  {  return *(this->_item);  }


private:

    bool            _isDone;
    size_t _currentBank;
    //vector<Iterator<Item>* > _refs;
    Iterator<Item>* _ref;
    size_t _nbBanks;
    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadProcessed;
    size_t _currentInternalBank;
	size_t _currentDataset;
	size_t _nbDatasets;

    Iterator<Item>* _mainref;
    void setMainref (Iterator<Item>* mainref)  { SP_SETATTR(mainref); }
};























struct SimkaSequenceFilter
{
	//u_int64_t _maxNbReads;
	//u_int64_t _maxNbReadsPerBank;
	//u_int64_t _nbReadProcessed;
	//CancellableIterator<Sequence>* _it;
	//int* _bankIndex;
	//int* _datasetIndex;


	SimkaSequenceFilter(size_t minReadSize, double minShannonIndex){
		//_maxNbReads = 0;
		//_nbReadProcessed = 0;
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

	//void setMaxReads(u_int64_t maxReads){
	//	_maxNbReads = maxReads;
	//}

	//void setIt(CancellableIterator<Sequence>* it){
	//	_it = it;
	//}

	bool operator() (Sequence& seq){

		//cout << seq.toString() << endl;
		//cout << _nbReadProcessed << endl;
		//if(_maxNbReads != 0){
		//	if(_nbReadProcessed >= _maxNbReads){
		//		_it->_cancel = true;
		//		return false;
		//	}
		//}

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
		//_nbReadProcessed += 1;

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



































/*********************************************************************
* ** SimkaCountProcessor
*********************************************************************/
template<size_t span>
class SimkaCountProcessorSimple{

private:

    size_t _nbBanks;
    size_t _nbNewBanks;
    size_t _bankOffset;
    size_t _kmerSize;
	//pair<CountNumber, CountNumber> _abundanceThreshold;
	//bool isAbundanceThreshold;

    SimkaStatistics* _stats;
    double _totalAbundance;

    u_int64_t _nbKmerCounted;
    double _minKmerShannonIndex;


    vector<u_int64_t> _sharedNewBanks;
    vector<u_int64_t> _sharedOldBanks;



public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

    SimkaCountProcessorSimple(SimkaStatistics* stats, size_t nbBanks, size_t nbNewBanks, size_t kmerSize, double minKmerShannonIndex) :
    _stats(stats)
    {

    	// We configure the vector for the N.(N+1)/2 possible pairs
    	//_countTotal.resize (_nbBanks*(_nbBanks+1)/2);

    	_nbBanks = nbBanks;
    	_nbNewBanks = nbNewBanks;
    	_kmerSize = kmerSize;
    	//_abundanceThreshold = abundanceThreshold;
    	_minKmerShannonIndex = minKmerShannonIndex;

    	//_localStats = new SimkaStatistics(_nbBanks, _stats._distanceParams);

    	_nbKmerCounted = 0;
    	//isAbundanceThreshold = _abundanceThreshold.first > 1 || _abundanceThreshold.second < 1000000;

    	_bankOffset = _nbBanks - _nbNewBanks;
    }


    void end(){
    }

    void process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& counts){
    	updateDistance(counts);
    }

    void updateDistance(const CountVector& counts){
    	_sharedNewBanks.clear();
    	_sharedOldBanks.clear();


		for(size_t i=0; i<_bankOffset; i++){
			if(counts[i]){
				_stats->_nbDistinctKmersPerDataset[i] += 1;
				_sharedOldBanks.push_back(i);
			}
		}
		for(size_t i=_bankOffset; i<counts.size(); i++){
			if(counts[i]){
				_stats->_nbDistinctKmersPerDataset[i] += 1;
				_sharedNewBanks.push_back(i);
			}
		}
		//for(size_t i=0; i<counts.size(); i++)
		//	if(counts[i]) _sharedBanks.push_back(i);




		if(_stats->_nbKmersPerDatasetPairs._matrix_squaredHalf.size() > 0){
			for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

				u_int64_t i = _sharedNewBanks[ii];
				//if(i == _nbBanks-1) continue;

				u_int64_t abundanceI = counts[i];

				int iMax = i-_bankOffset-1;
				if(iMax >= 0){
					for(int i2=0; i2<=iMax; i2++){
						u_int64_t abundanceJ = counts[_bankOffset+i2];

						if(_stats->_nbDistinctKmersPerDataset[i] + _stats->_nbDistinctKmersPerDataset[_bankOffset+i2] - _stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[i2][iMax-i2] >= _stats->_sketchSize){ //_minHashUnion[i] is _skecthSize
							continue;
						}

						_stats->_nbKmersPerDatasetPairs._matrix_squaredHalf[i2][iMax-i2] += abundanceI + abundanceJ;
					}
				}

				int iLocal = i-_bankOffset;
				if(iLocal < _stats->_nbKmersPerDatasetPairs._matrix_squaredHalf.size()){
					for(int j=0; j<_stats->_nbKmersPerDatasetPairs._matrix_squaredHalf[iLocal].size(); j++){

						size_t jGlobal = iLocal+j+1+_bankOffset;

						if(_stats->_nbDistinctKmersPerDataset[i] + _stats->_nbDistinctKmersPerDataset[jGlobal] - _stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[iLocal][j] >= _stats->_sketchSize){ //_minHashUnion[i] is _skecthSize
							continue;
						}

						u_int64_t abundanceJ = counts[jGlobal];
						_stats->_nbKmersPerDatasetPairs._matrix_squaredHalf[iLocal][j] += abundanceI + abundanceJ;
					}
				}

			}
		}

		for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

			u_int64_t i = _sharedNewBanks[ii];
			u_int64_t abundanceI = counts[i];

			for(size_t jj=0; jj<_bankOffset; jj++){

				u_int64_t j = _sharedOldBanks[jj];
				u_int64_t abundanceJ = counts[j];

				if(_stats->_nbDistinctKmersPerDataset[i] + _stats->_nbDistinctKmersPerDataset[j] - _stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[i-_bankOffset][j] >= _stats->_sketchSize){ //_minHashUnion[i] is _skecthSize
					continue;
				}

				_stats->_nbKmersPerDatasetPairs._matrix_squaredHalf[i-_bankOffset][j] += abundanceI + abundanceJ;

			}
		}


		/*
		for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

			u_int64_t i = _sharedNewBanks[ii];
			u_int64_t abundanceI = counts[i];

			_minHashUnion[i] += 1;

			for(size_t jj=ii+1; jj<_sharedNewBanks.size(); jj++){

				u_int64_t j = _sharedNewBanks[jj] - _bankOffset;
				u_int64_t abundanceJ = counts[j];

				_minHashUnionCross[i][j] += 1;
			}
		}*/

		updateDistanceDefault(counts);
    }

	void updateDistanceDefault(const CountVector& counts){

		/*
		cout << _stats->_brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		for(size_t i=0; i<_stats->_brayCurtisNumerator._matrix_squaredHalf.size(); i++){
			cout << i << endl;
			for(size_t j=0; j<_stats->_brayCurtisNumerator._matrix_squaredHalf[i].size(); j++){
				cout << "\t" << j << endl;
			}
		}*/

		if(_stats->_brayCurtisNumerator._matrix_squaredHalf.size() > 0){
			for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

				u_int64_t i = _sharedNewBanks[ii];
				if(i == _nbBanks-1) continue;

				u_int64_t abundanceI = counts[i];
				u_int64_t jOffset = _stats->_brayCurtisNumerator._matrix_squaredHalf.size() - _stats->_brayCurtisNumerator._matrix_squaredHalf[i-_bankOffset].size();// + 1;

				for(size_t jj=ii+1; jj<_sharedNewBanks.size(); jj++){

					u_int64_t j = _sharedNewBanks[jj];

					//cout << _stats->_nbSolidDistinctKmersPerBank[i] << endl;
					if(_stats->_nbDistinctKmersPerDataset[i] + _stats->_nbDistinctKmersPerDataset[j] - _stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[i-_bankOffset][j-_bankOffset-jOffset-1] >= _stats->_sketchSize){ //_minHashUnion[i] is _skecthSize
						continue;
					}
					//cout << _minHashUnion[i] << " " <<  _minHashUnion[j] << " " << _minHashUnionCross[i][j] << "    " << (_minHashUnion[i] + _minHashUnion[j] - _minHashUnionCross[i][j]) << "  " << _stats->_nbSolidDistinctKmersPerBank[i] << endl;

					//size_t symetricIndex = j + ((_nbBanks-1)*i) - (i*(i-1)/2);

					u_int64_t abundanceJ = counts[j];

					//_stats->_matrixNbSharedKmers[i][j] += counts[i];
					//_stats->_matrixNbSharedKmers[j][i] += counts[j];
					//_stats->_matrixNbDistinctSharedKmers[i][j] += 1;


					//cout << "\t    " << _stats->_brayCurtisNumerator._matrix_squaredHalf.size() << " " << i-_bankOffset << endl;
					//cout << "\t        " << _stats->_brayCurtisNumerator._matrix_squaredHalf[i-_bankOffset].size() << " " << (j-_bankOffset-jOffset-1) << endl;
					//cout << i << " " << j << "     " << _stats->_brayCurtisNumerator._matrix_squaredHalf.size() << " " << _stats->_brayCurtisNumerator._matrix_squaredHalf[i].size() << endl;
					//cout << i-jOffset << " " << j-jOffset-1 << endl;
					//cout << jOffset << endl;
					_stats->_brayCurtisNumerator._matrix_squaredHalf[i-_bankOffset][j-_bankOffset-jOffset-1] += min(abundanceI, abundanceJ);
					_stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[i-_bankOffset][j-_bankOffset-jOffset-1] += 1;

				}
			}
		}

		/*
		cout << "---" << endl;
		for(size_t i=0; i<counts.size(); i++){
			cout << counts[i] << "\t";
		}
		cout << endl;
		for(size_t i=0; i<_stats->_matrixNbDistinctSharedKmers.size(); i++){
			for(size_t j=0; j<_stats->_matrixNbDistinctSharedKmers[i].size(); j++){
				cout << _stats->_matrixNbDistinctSharedKmers[i][j] << "\t";
			}
			cout << endl;
		}
		cout << "---" << endl;*/

		//cout << "---" << endl;

		for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

			u_int64_t i = _sharedNewBanks[ii];
			u_int64_t abundanceI = counts[i];

			for(size_t jj=0; jj<_sharedOldBanks.size(); jj++){

				u_int64_t j = _sharedOldBanks[jj];

				if(_stats->_nbDistinctKmersPerDataset[i] + _stats->_nbDistinctKmersPerDataset[j] - _stats->_matrixNbDistinctSharedKmers._matrix_squaredHalf[i-_bankOffset][j] >= _stats->_sketchSize){
					continue;
				}

				u_int64_t abundanceJ = counts[j];

				_stats->_brayCurtisNumerator._matrix_rectangular[i-_bankOffset][j] += min(abundanceI, abundanceJ);
				_stats->_matrixNbDistinctSharedKmers._matrix_rectangular[i-_bankOffset][j] += 1;

			}
		}

	}


};




















class SimkaCounterBuilderMerge
{
public:

    /** Constructor.
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
	SimkaCounterBuilderMerge (CountVector& abundancePerBank)  :  _abundancePerBank(abundancePerBank)  {}

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank, CountNumber abundance)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= abundance;
    	//cout << _abundancePerBank.size() << " " << idxBank << " " << abundance << endl;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank, CountNumber abundance)  {
    	//cout << _abundancePerBank.size() << " " << idxBank << " " << abundance << endl;
    	_abundancePerBank [idxBank] += abundance;
    }

    /** Set the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    //void set (CountNumber val, size_t idxBank=0)  {  _abundancePerBank [idxBank] = val;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    //CountNumber operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** */
    //const CountVector& get () const { return _abundancePerBank; }

    void print(const string& kmer){
		cout << kmer << ": ";
    	for(size_t i=0; i<size(); i++){
    		cout << _abundancePerBank[i] << " ";
    	}
    	cout << endl;
    }

    CountVector& _abundancePerBank;
};



#endif /* GATB_SIMKA_SRC_CORE_UTILS_HPP_ */
