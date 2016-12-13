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

    //vector<size_t> _banksOks;

    vector<u_int64_t> _sharedNewBanks;
    vector<u_int64_t> _sharedOldBanks;

	typedef std::pair<double, CountVector> chi2val_Abundances;
	struct _chi2ValueSorterFunction { bool operator() (chi2val_Abundances l,chi2val_Abundances r) { return r.first < l.first; } } ;
	std::priority_queue< chi2val_Abundances, vector<chi2val_Abundances>, _chi2ValueSorterFunction> _chi2ValueSorter;
	size_t _maxChi2Values;

public:

    typedef typename Kmer<span>::Type  Type;
    //typedef typename Kmer<span>::Count Count;

    SimkaCountProcessorSimple(SimkaStatistics* stats, size_t nbBanks, size_t nbNewBanks, size_t kmerSize, double minKmerShannonIndex) :
    _stats(stats)
    {

    	_maxChi2Values = 1000;
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
		#ifdef CHI2_TEST

			size_t nbValues = _chi2ValueSorter.size();
			for(size_t i=0; i<nbValues; i++){
				double val = _chi2ValueSorter.top().first;
				CountVector counts = _chi2ValueSorter.top().second;


				//cout << val << endl;
				//for(size_t j=0; j<counts.size(); j++){
				//	cout << counts[j] << " ";
				//}
				//cout << endl;

				updateDistance(counts);

				_chi2ValueSorter.pop();
			}
		#endif
    }

    void process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& counts){

    	//cout << kmer.toString(_kmerSize) << endl;
    	//for(size_t i=0; i<counts.size(); i++){
    	//	cout << counts[i] << " ";
    	//}
    	//cout << endl;



    	/* A DEPLACER PENDANT LE COMPTAGE DES KMERS
    	if(_minKmerShannonIndex != 0){
    		double shannonIndex = getShannonIndex(kmer);
    		if(shannonIndex < _minKmerShannonIndex){
    			return;
    		}
    	}*/

#ifdef CHI2_TEST
    	float X2j = 0;

    	_totalAbundance = 0;
    	for(size_t i=0; i<counts.size(); i++){
    		_totalAbundance += counts[i];
    	}

    	for(size_t i=0; i<counts.size(); i++){
    		X2j += pow((counts[i]/_totalAbundance - _stats->_datasetNbReads[i]/_stats->_totalReads), 2) / (_stats->_datasetNbReads[i] / (_stats->_totalReads*_totalAbundance));
    	}

    	//std::chi_squared_distribution<double> distribution(_nbBanks-1);
    	//double pvalue = chisqr(_nbBanks-1, X2j);

    	/*
    	if(lala> 100){
        	for(size_t i=0; i<_chi2ValueSorter.size(); i++){
        		double val = _chi2ValueSorter.top();
        		_chi2ValueSorter.pop();
        		cout << val << endl;
        	}

    		return;
    	}*/

    	//cout << X2j << endl;


    	if(_chi2ValueSorter.size() > _maxChi2Values){

        	if(X2j > _chi2ValueSorter.top().first){
            	_chi2ValueSorter.push(pair<double, CountVector>(X2j, counts));
        		_chi2ValueSorter.pop();
        	}

    	}
    	else{
        	_chi2ValueSorter.push(pair<double, CountVector>(X2j, counts));
    	}



    	//cout << _chi2ValueSorter.size() << "    " << X2j << "   " << _chi2ValueSorter.top() << endl;
    	//cout <<  X2j << "    " << pvalue << endl;

    	return;
    	/*
    	cout << kmer.toString(_kmerSize) << "  [";
    	for(size_t i=0; i<counts.size(); i++)
    			cout << counts[i] << " ";
    	cout << "]    ";
    	cout <<  X2j << "    " << pvalue << endl;
*/

    	//cout
    	//cout << X2j << endl;
    	//if(_totalAbundance == 1){
    	//	cout << X2j << endl;
    	//}
    	//if(pvalue > 0.01) return;
#endif
    	/*
    	//for(size_t i=0; i<_datasetNbReads.size(); i++)
    	//	cout << i << " " << _datasetNbReads[i] << endl;

    	//cout << _totalReads << " " << _totalAbundance << endl;
    	//float Ri = 500000;
    	//float Rtotal = Ri * _nbBanks;
    	//float Ntotal = _totalAbundance;
    	float X2j = 0;
    	for(size_t i=0; i<counts.size(); i++){

    		float Ni = counts[i];

    		X2j += pow((Ni/_totalAbundance - _datasetNbReads[i]/_totalReads), 2) / (_datasetNbReads[i] / (_totalReads*_totalAbundance));
    	}

    	//cout << X2j << endl;
    	//if(_totalAbundance == 1){
    	//	cout << X2j << endl;
    	//}
    	if(X2j <= (_nbBanks-1)*1.7) return false;
    	*/

    	//cout << X2j << endl;


    	//cout << kmer.toString(31) << endl;

    	//cout << endl;

    	//if(_progress){ //Simka_min
    	//	_stats->_nbSolidKmers += 1;
    	//	computeStats(counts);
    	//}
    	//else{

    	updateDistance(counts);

    	//else
    	//	computeStats(counts);


    	//_stats->_nbSolidKmers += 1;
    }

    void updateDistance(const CountVector& counts){
    	_sharedNewBanks.clear();
    	_sharedOldBanks.clear();


		for(size_t i=0; i<_bankOffset; i++){
			if(counts[i]) _sharedOldBanks.push_back(i);
		}
		for(size_t i=_bankOffset; i<counts.size(); i++){
			if(counts[i]) _sharedNewBanks.push_back(i);
		}
		//for(size_t i=0; i<counts.size(); i++)
		//	if(counts[i]) _sharedBanks.push_back(i);

		updateDistanceDefault(counts);

    	if(_stats->_computeSimpleDistances)
    		updateDistanceSimple(counts);

    	if(_stats->_computeComplexDistances)
    		updateDistanceComplex(counts);
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

		for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

			u_int64_t i = _sharedNewBanks[ii];
			u_int64_t abundanceI = counts[i];
			u_int64_t jOffset = _stats->_brayCurtisNumerator._matrix_squaredHalf.size() - _stats->_brayCurtisNumerator._matrix_squaredHalf[i-_bankOffset].size();// + 1;

			for(size_t jj=ii+1; jj<_sharedNewBanks.size(); jj++){

				u_int64_t j = _sharedNewBanks[jj];
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

				//u_int64_t i = _sharedNewBanks[ii];
				u_int64_t j = _sharedOldBanks[jj];
				//size_t symetricIndex = j + ((_nbBanks-1)*i) - (i*(i-1)/2);

				//u_int64_t abundanceI = counts[i];
				u_int64_t abundanceJ = counts[j];

				//_stats->_matrixNbSharedKmers[i][j] += counts[i];
				//_stats->_matrixNbSharedKmers[j][i] += counts[j];
				//_stats->_matrixNbDistinctSharedKmers[i][j - _bankOffset] += 1;

				//cout << i << " " << j << "    " << (j + ((_nbBanks-1)*i) - (i*(i-1)/2)) << endl;
				//cout << _stats->_brayCurtisNumerator._matrix_rectangular.size() << "      " << i << " " << j << endl;
				//cout << "\t" << _stats->_brayCurtisNumerator._matrix_rectangular.size() << " " << i-_bankOffset << endl;
				//cout << "\t" << _stats->_brayCurtisNumerator._matrix_rectangular[i-_bankOffset].size() << " " << j << endl;
				_stats->_brayCurtisNumerator._matrix_rectangular[i-_bankOffset][j] += min(abundanceI, abundanceJ);
				_stats->_matrixNbDistinctSharedKmers._matrix_rectangular[i-_bankOffset][j] += 1;

				//(counts.size()-1) - (j - _bankOffset) - _bankOffset
				//cout << _stats->_brayCurtisNumerator.size() << " " << _stats->_brayCurtisNumerator[0].size() << "     " << i << " " << j << endl;
				//cout << i << " " << j << endl;
			}
		}

	}


	void updateDistanceSimple(const CountVector& counts){


		for(size_t ii=0; ii<_sharedNewBanks.size(); ii++){

			u_int64_t i = _sharedNewBanks[ii];
			u_int64_t abundanceI = counts[i];

			for(size_t jj=ii+1; jj<_sharedNewBanks.size(); jj++){

				u_int64_t j = _sharedNewBanks[jj] - _bankOffset;
				u_int64_t abundanceJ = counts[j];

				//cout << _stats->_chord_sqrt_N2[i] << endl;
				//_stats->_chord_NiNj[i][j] += abundanceI * abundanceJ;
				_stats->_chord_NiNj[i][j] += abundanceI * abundanceJ;
				_stats->_hellinger_SqrtNiNj[i][j] += sqrt(abundanceI * abundanceJ);
				_stats->_kulczynski_minNiNj[i][j] += min(abundanceI, abundanceJ);
			}
		}


		for(size_t ii=0; ii<_sharedOldBanks.size(); ii++){

			u_int64_t i = _sharedOldBanks[ii];
			u_int64_t abundanceI = counts[i];

			for(size_t jj=0; jj<_sharedNewBanks.size(); jj++){

				u_int64_t j = _sharedNewBanks[jj] - _bankOffset;
				u_int64_t abundanceJ = counts[j];

				_stats->_chord_NiNj[i][j] += abundanceI * abundanceJ;
				_stats->_hellinger_SqrtNiNj[i][j] += sqrt(abundanceI * abundanceJ);
				_stats->_kulczynski_minNiNj[i][j] += min(abundanceI, abundanceJ);
			}
		}


	}

	void updateDistanceComplex(const CountVector& counts){


		//_sharedBanks.clear();

		//for(size_t i=0; i<counts.size(); i++)
		//	if(counts[i]) _sharedBanks.push_back(i);

		double xi = 0;
		double xj = 0;
		double d1 = 0;
		double d2 = 0;

		for(size_t i=0; i<counts.size(); i++){
			if(counts[i]){
				for(size_t j=i+1; j<counts.size(); j++){

					//In this loop we know that (abundanceI > 0)
					double abundanceI = counts[i];
					double abundanceJ = counts[j];

					if(abundanceJ){


						//_stats->_matrixNbSharedKmers[i][j] += abundanceI;
						//_stats->_matrixNbSharedKmers[j][i] += abundanceJ;
						//_stats->_matrixNbDistinctSharedKmers[i][j] += 1;
						//_stats->_chord_NiNj[i][j] += abundanceI * abundanceJ;
						//_stats->_chord_NiNj[i][j] += (abundanceI * abundanceJ) / (_stats->_chord_sqrt_N2[i]*_stats->_chord_sqrt_N2[j]);
						//_stats->_hellinger_SqrtNiNj[i][j] += sqrt(abundanceI * abundanceJ);
						//_stats->_kulczynski_minNiNj[i][j] += min(abundanceI, abundanceJ);


						double yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
						double xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
						xi = (double)abundanceI / _stats->_nbSolidKmersPerBank[i];
						d1 = xi * log((2*xY) / (xY + yX));


						//xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
						//yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
						xj = (double)abundanceJ / _stats->_nbSolidKmersPerBank[j];
						d2 = xj * log((2*yX) / (xY + yX));
					}
					else{
						d2 = 0;
						double yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
						double xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
						xi = (double)abundanceI / _stats->_nbSolidKmersPerBank[i];
						d1 = xi * log((2*xY) / (xY + yX));
					}

					/*
					if(abundanceI){
						double yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
						double xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
						xi = (double)abundanceI / _stats->_nbSolidKmersPerBank[i];
						d1 = xi * log((2*xY) / (xY + yX));
					}
					else{
						d1 = 0;
					}

					if(abundanceJ){
						double xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
						double yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
						xj = (double)abundanceJ / _stats->_nbSolidKmersPerBank[j];
						d2 = xj * log((2*yX) / (xY + yX));
					}
					else{
						d2 = 0;
					}*/

					_stats->_kullbackLeibler[i][j] += d1 + d2;

					_stats->_canberra[i][j] += abs(abundanceI - abundanceJ) / (abundanceI + abundanceJ);
					//_stats->_brayCurtisNumerator[i][j] += abs(abundanceI - abundanceJ);
					_stats->_whittaker_minNiNj[i][j] += abs((int)((u_int64_t)(abundanceI*_stats->_nbSolidKmersPerBank[j]) - (u_int64_t)(abundanceJ*_stats->_nbSolidKmersPerBank[i])));

					//cout << _stats->_nbSolidKmersPerBank[i] << endl;


				}
			}
			else{

				//Here, we know that (abundanceI == 0)

				for(size_t jj=0; jj<_sharedNewBanks.size(); jj++){

					u_int16_t j = _sharedNewBanks[jj];
					if(i > j) continue;

					double abundanceI = counts[i];
					double abundanceJ = counts[j];

					d1 = 0;
					double xY = abundanceI * _stats->_nbSolidKmersPerBank[j];
					double yX = abundanceJ * _stats->_nbSolidKmersPerBank[i];
					xj = (double)abundanceJ / _stats->_nbSolidKmersPerBank[j];
					d2 = xj * log((2*yX) / (xY + yX));

					_stats->_kullbackLeibler[i][j] += d1 + d2;

					_stats->_canberra[i][j] += abs(abundanceI - abundanceJ) / (abundanceI + abundanceJ);
					//_stats->_brayCurtisNumerator[i][j] += abs(abundanceI - abundanceJ);
					//cout << _stats->_nbSolidKmersPerBank[i] << endl;

					_stats->_whittaker_minNiNj[i][j] += abs((int)((u_int64_t)(abundanceI*_stats->_nbSolidKmersPerBank[j]) - (u_int64_t)(abundanceJ*_stats->_nbSolidKmersPerBank[i])));

				}
			}
		}


	}

	//inline bool isSolidVector(const CountVector& counts);
	double getShannonIndex(const Type&  kmer){
		float index = 0;
		//float freq [5];

		vector<float> _freqs(4, 0);

		//char* seqStr = seq.getDataBuffer();

	    for (size_t i=0; i<_kmerSize; i++){
	    	_freqs[kmer[i]] += 1.0;
	    	//seq[sizeKmer-i-1] = bin2NT [(*this)[i]];
	    }

		// Frequency of each letter (A, C, G, T or N)
		//for(size_t i=0; i < seq.size(); i++)
		//	_freqs[nt2binTab[(unsigned char)seq[i]]] += 1.0;

		// Shannon index calculation
		for (size_t i=0; i<_freqs.size(); i++){
			_freqs[i] /= (float) _kmerSize;
			if (_freqs[i] != 0)
				index += _freqs[i] * log (_freqs[i]) / log(2);
		}
		return abs(index);
	}

	double approx_gamma(double Z)
	{
	    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
	    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

	    double D = 1.0 / (10.0 * Z);
	    D = 1.0 / ((12 * Z) - D);
	    D = (D + Z) * RECIP_E;
	    D = pow(D, Z);
	    D *= sqrt(TWOPI / Z);

	    return D;
	}

	static double igf(double S, double Z)
	{
	    if(Z < 0.0)
	    {
		return 0.0;
	    }
	    double Sc = (1.0 / S);
	    Sc *= pow(Z, S);
	    Sc *= exp(-Z);

	    double Sum = 1.0;
	    double Nom = 1.0;
	    double Denom = 1.0;

	    for(int I = 0; I < 200; I++)
	    {
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	    }

	    return Sum * Sc;
	}

	double chisqr(int Dof, double Cv)
	{
	    if(Cv < 0 || Dof < 1)
	    {
	        return 0.0;
	    }
	    double K = ((double)Dof) * 0.5;
	    double X = Cv * 0.5;
	    if(Dof == 2)
	    {
		return exp(-1.0 * X);
	    }

	    double PValue = igf(K, X);
	    //if(isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
	    //{
	    //    return 1e-14;
	    //}

	    PValue /= approx_gamma(K);
	    //PValue /= tgamma(K);

	    return PValue;
	    //return (1.0 - PValue);
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
