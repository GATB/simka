/*
 * SimkaCommons.h
 *
 *  Created on: 24 juin 2017
 *      Author: gbenoit
 */

#ifndef SIMKA1_4_SRC_CORE_SIMKACOMMONS_HPP_
#define SIMKA1_4_SRC_CORE_SIMKACOMMONS_HPP_

#include <thread>

const string STR_SIMKA_SOLIDITY_PER_DATASET = "-solidity-single";
const string STR_SIMKA_MAX_READS = "-max-reads";
const string STR_SIMKA_MIN_READ_SIZE = "-min-read-size";
const string STR_SIMKA_MIN_READ_SHANNON_INDEX = "-min-shannon-index";
const string STR_SIMKA_MIN_KMER_SHANNON_INDEX = "-kmer-shannon-index";
const string STR_KMER_PER_READ = "-kmer-per-read";
const string STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES= "-simple-dist";
const string STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES = "-complex-dist";
const string STR_SIMKA_KEEP_TMP_FILES = "-keep-tmp";
const string STR_SIMKA_COMPUTE_DATA_INFO = "-data-info";



class SimkaCommons {
public:
	SimkaCommons();
	virtual ~SimkaCommons();


	static void checkInputValidity(const string& outputDirTemp, const string& inputFilename, u_int64_t& nbDatasets){

		if(!System::file().doesExist(inputFilename)){
			cout << "ERROR: Input does not exists (" + inputFilename + ")" << endl;
			exit(1);
		}

		nbDatasets = 0;
		bool error = false;

		//string inputDir = _outputDirTemp; // + "/input/";
		ifstream inputFile(inputFilename.c_str());

		//ofstream outputFileIds(_outputFilenameIds.c_str(), ios::binary);
		//_banksInputFilename =  inputDir + "__input_simka__"; //_inputFilename + "_dsk_dataset_temp__";
		//IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		string line;
		string linePart;
		vector<string> lineIdDatasets;
		vector<string> linepartPairedDatasets;
		vector<string> linepartDatasets;

		//string bankFileContents = "";

		u_int64_t lineIndex = 0;
		u_int64_t bankIdBytePos = 0;

		while(getline(inputFile, line)){

			line.erase(std::remove(line.begin(),line.end(),' '),line.end());
			if(line == "") continue;

			//cout << line << endl;
			lineIdDatasets.clear();
			linepartPairedDatasets.clear();
			//vector<string> filenames;

			stringstream lineStream(line);
			while(getline(lineStream, linePart, ':')){
				lineIdDatasets.push_back(linePart);
			}

			string bankId = lineIdDatasets[0];
			string linePairedDatasets = lineIdDatasets[1];

			stringstream linePairedDatasetsStream(linePairedDatasets);
			while(getline(linePairedDatasetsStream, linePart, ';')){
				linepartPairedDatasets.push_back(linePart);
			}

			string subBankFilename = outputDirTemp + bankId;
			IFile* subBankFile = System::file().newFile(subBankFilename, "wb");
			//cout << subBankFile->getPath() << endl;
			string subBankContents = "";
			//_nbBankPerDataset.push_back(linepartPairedDatasets.size());

			for(size_t i=0; i<linepartPairedDatasets.size(); i++){
				string lineDatasets = linepartPairedDatasets[i];

				linepartDatasets.clear();

				stringstream lineDatasetsStream(lineDatasets);
				while(getline(lineDatasetsStream, linePart, ',')){
					linepartDatasets.push_back(linePart);
					//cout << "\t" << linePart << endl;
				}

				//bankFileContents += linepartDatasets[0] + "\n";


				for(size_t i=0; i<linepartDatasets.size(); i++){
					string filename = linepartDatasets[i];
					if(filename.at(0) == '/'){
						subBankContents +=  filename + "\n";
					}
					else{
						string dir = System::file().getRealPath(inputFilename);
						dir = System::file().getDirectory(dir);
						subBankContents +=  dir + "/" + filename + "\n";
					}
				}

			}

			subBankContents.erase(subBankContents.size()-1);
			subBankFile->fwrite(subBankContents.c_str(), subBankContents.size(), 1);
			subBankFile->flush();
			delete subBankFile;

			//bankFileContents += inputDir + "/" + bankId + "\n";
			lineIndex += 1;


			try{
				IBank* bank = Bank::open(subBankFilename);
				LOCAL(bank);
				nbDatasets += 1;
			}
			catch (Exception& e){
				cerr << "ERROR: Can't open dataset: " << bankId << endl;
				error = true;
			}

			System::file().remove(subBankFilename);

		}


		inputFile.close();

		if(error) exit(1);

	}


};










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






template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:


	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads, size_t nbDatasets) : BankDelegate (ref), _ref2(0), _filter(filter)  {
		_maxReads = maxReads;
		_nbDatasets = nbDatasets;
		setRef2(_ref->iterator ());
	}



	~SimkaPotaraBankFiltered(){

	    std::vector<Iterator<Sequence>*> itBanks =  _ref2->getComposition();
	    for(size_t i=0; i<itBanks.size(); i++){
	    	delete itBanks[i];
	    }

	    //_ref2->
		setRef2(0);
	}

    Iterator<Sequence>* iterator ()
    {
        return new SimkaInputIterator<Sequence, Filter> (_ref2, _nbDatasets, _maxReads, _filter);

    }

private:


    Iterator<Sequence>* _ref2;
    void setRef2 (Iterator<Sequence>* ref2)  { SP_SETATTR(ref2); }

    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadToProcess;
    size_t _datasetId;
    size_t _nbDatasets;
};




















#endif /* SIMKA1_4_SRC_CORE_SIMKACOMMONS_H_ */
