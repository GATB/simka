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

#include "SimkaAlgorithm.hpp"
//#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;









/********************************************************************************/
/** \brief Iterator that can be cancelled at some point during iteration
 *
 * This iterator iterates a referred iterator and will finish:
 *      - when the referred iterator is over
 *   or - when the cancel member variable is set to true
 *
 */
template <class Item, typename Filter> class SimkaInputIterator : public Iterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] initRef : will call 'first' on the reference if true
     */
	SimkaInputIterator(vector<Iterator<Item>* >& refs, size_t nbBanks, u_int64_t maxReads, Filter filter)
        :  _refs(refs), _ref(refs[0]), _filter(filter) {
		_isDone = true;
		_nbBanks = nbBanks;
		cout << _nbBanks << endl;
		_maxReads = maxReads;
		_nbReadProcessed = 0;
		_currentBank = 0;
		_currentInternalBank = 0;
	}

    /** \copydoc  Iterator::first */
    void first()
    {


        _ref->first();


        _isDone = _ref->isDone();


        while (!_ref->isDone() && _filter(_ref->item())==false)
        	_ref->next();


    	*(this->_item) = _ref->item();

    	//_ref->next();
        /** IMPORTANT : we need to copy the referred item => we just can't rely on a simple
         * pointer (in case of usage of Dispatcher for instance where a buffer of items is kept
         * and could be wrong if the referred items doesn't exist any more when accessing the
         * buffer).
         * TODO doc: I get it, but why is it done only in this iterator and not other iterators like FilterIterator?*/
        //if (!_isDone)  { *(this->_item) = _ref->item(); }

        //cout << "HAA   " << _currentBank << endl;
    }

    bool isFinished(){
    	if(_currentBank == _refs.size()-1){
    		_isDone = true;
    		return true;
    	}
    	return false;
    }

    void nextDataset(){
    	//cout << "next dataset "<< endl;
    	while(_currentInternalBank < _nbBanks){
        	_currentBank += 1;
    		_currentInternalBank += 1;
    	}
    	_currentInternalBank = 0;
    	_nbReadProcessed = 0;

		if(isFinished()){
			return;
		}

    	nextBank();
    }

    void nextBank(){
    	//cout << "next bank "<< endl;
    	_currentInternalBank += 1;
    	if(_currentInternalBank == _nbBanks){
    		nextDataset();
    	}
    	else{
        	_isDone = false;
        	_currentBank += 1;
        	_ref = _refs[_currentBank];
        	first();
    	}

    }





    /** \copydoc  Iterator::next */
    void next()
    {

    	_ref->next();

        _isDone = _ref->isDone();

        while (!_ref->isDone() && _filter(_ref->item())==false)
        	_ref->next();

        //_isDone = _ref->isDone();

        //if (!_isDone)  {
        	*(this->_item) = _ref->item();
        	_nbReadProcessed += 1;

        	//cout << &_ref->item() << endl;
        //}

    	//cout << _nbReadProcessed << "  " << _maxReads << "    " << _refs.size() << endl;


        if(_isDone){
    		if(isFinished())
    			return;
    		else
    			nextBank();

        }
        else{
        	//*(this->_item) = _ref->item();
        }

    	if(_nbReadProcessed >= _maxReads){
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
    vector<Iterator<Item>* > _refs;
    Iterator<Item>* _ref;
    size_t _nbBanks;
    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadProcessed;
    size_t _currentInternalBank;
};



template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:

	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads, size_t nbDatasets) : BankDelegate (ref), _filter(filter)  {
		//_nbReadsPerDataset = nbReadsPerDataset;
		_maxReads = maxReads;
		_nbDatasets = nbDatasets;
	}


    Iterator<Sequence>* iterator ()
    {

        Iterator<Sequence>* it = _ref->iterator ();
        std::vector<Iterator<Sequence>*> iterators = it->getComposition();
        return new SimkaInputIterator<Sequence, Filter> (iterators, _nbDatasets, _maxReads, _filter);
    	//return filterIt;

    }

private:

	//vector<u_int64_t> _nbReadsPerDataset;
    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadToProcess;
    size_t _datasetId;
    size_t _nbDatasets;
};


class SimkaCount : public Tool
{
public:

	SimkaCount () : Tool ("SimkaCount")
    {
        //getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "output file",           true));
        //getParser()->push_back (new OptionOneParam (STR_ID,   "dataset id", true));
        //getParser()->push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size", true));
        getParser()->push_back (new OptionOneParam ("-out-tmp-simka",   "tmp output", true));
        getParser()->push_back (new OptionOneParam ("-bank-name",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MAX_READS,   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-datasets",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));

        getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);
        if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    }

    void execute ()
    {


    	//size_t datasetId =  getInput()->getInt(STR_ID);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	//cout << kmerSize << endl;

    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	string bankName =  getInput()->getStr("-bank-name");
    	size_t minReadSize =  getInput()->getInt(STR_SIMKA_MIN_READ_SIZE);
    	double minReadShannonIndex =  getInput()->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
    	u_int64_t maxReads =  getInput()->getInt(STR_SIMKA_MAX_READS);
    	size_t nbDatasets =   getInput()->getInt("-nb-datasets");

    	Parameter params(*this, kmerSize, outputDir, bankName, minReadSize, minReadShannonIndex, maxReads, nbDatasets);

        Integer::apply<Functor,Parameter> (kmerSize, params);



		//SimkaBankId* bank = new SimkaBankId(_banks, i);
		//cout << config._nb_partitions << endl;
		//KmerCountCompressor<span>* kmerCountCompressor = new KmerCountCompressor<span>(outputDir, config._nb_partitions, 1);

		//SimkaCompProcessor<span>* processor = new SimkaCompProcessor<span>(kmerCountCompressor);
		//vector<ICountProcessor<span>*> procs;
		//procs.push_back(processor);

		//algo.addProcessor(processor);

		//algo.execute();

		//delete kmerCountCompressor;
		//itBanks[i]->


        // We get a handle on the HDF5 storage object.
        // Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
        //Storage* storage = StorageFactory(DSK::getStorageMode()).load (getInput()->getStr(STR_URI_FILE));
        //LOCAL (storage);

        //string kmerSizeStr = storage->getGroup("params").getProperty ("kmer_size");

        //if (kmerSizeStr.empty())  { throw Exception ("unable to get the kmer size"); }

        //size_t kmerSize = atoi (kmerSizeStr.c_str());

    }


    struct Parameter
    {
        Parameter (SimkaCount& tool, size_t kmerSize, string outputDir, string bankName, size_t minReadSize, double minReadShannonIndex, u_int64_t maxReads, size_t nbDatasets) :
        	tool(tool), kmerSize(kmerSize), outputDir(outputDir), bankName(bankName), minReadSize(minReadSize), minReadShannonIndex(minReadShannonIndex), maxReads(maxReads), nbDatasets(nbDatasets)  {}
        SimkaCount& tool;
        //size_t datasetId;
        size_t kmerSize;
        string outputDir;
        string bankName;
        size_t minReadSize;
        double minReadShannonIndex;
        u_int64_t maxReads;
        size_t nbDatasets;
    };

    template<size_t span> struct Functor  {

    	void operator ()  (Parameter p){


			IProperties* props = p.tool.getInput();


			Configuration* config = new Configuration();
			{
				Repartitor* repartitor = new Repartitor();
				LOCAL(repartitor);

				{
					Storage* storage = StorageFactory(STORAGE_HDF5).load (p.outputDir + "/" + "config.h5");
					LOCAL (storage);
					config->load(storage->getGroup(""));
					repartitor->load(storage->getGroup(""));
				}

				config->_abundanceUserNb = 1;
				config->_abundance.clear();
				CountRange range(props->getInt(STR_KMER_ABUNDANCE_MIN), 100000);
				config->_abundance.push_back(range);

				//delete storage;
				/*
				config._kmerSize = p.kmerSize;
				config._minim_size = 8;
				config._max_disk_space = 0;
				config._max_memory = props->getInt(STR_MAX_MEMORY);
				config._nbCores = props->getInt(STR_NB_CORES);
				config._nb_partitions_in_parallel = config._nbCores;*/
				/*
				size_t      _kmerSize;
				size_t      _minim_size;
				size_t      _repartitionType;
				size_t      _minimizerType;

				tools::misc::KmerSolidityKind _solidityKind;

				u_int64_t   ;
				u_int32_t   _max_memory;

				size_t      _nbCores;
				size_t      _nb_partitions_in_parallel;
				size_t      _partitionType;

				std::vector<tools::misc::CountRange>  _abundance;
				size_t _abundanceUserNb;*/

				string tempDir = p.outputDir + "/temp/" + p.bankName;
				System::file().mkdir(tempDir, -1);
				//cout << i << endl;
				//string outputDir = p.outputDir + "/comp_part" + to_string(p.datasetId) + "/";

				//cout << "\tinput: " << p.outputDir + "/input/" + p.bankName << endl;
				IBank* bank = Bank::open(p.outputDir + "/input/" + p.bankName);

				SimkaSequenceFilter sequenceFilter(p.minReadSize, p.minReadShannonIndex);
				IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, p.maxReads, p.nbDatasets);
				// = new SimkaPotaraBankFiltered(bank)
				LOCAL(filteredBank);
				LOCAL(bank);

				Storage* solidStorage = 0;

				string solidsName = p.outputDir + "/solid/" +  p.bankName + ".h5";
				bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
				solidStorage = StorageFactory(STORAGE_HDF5).create (solidsName, true, autoDelete);
				LOCAL(solidStorage);


				//props->add(1, STR_HISTOGRAM_MAX, "0");
				//props->add(1, STR_KMER_ABUNDANCE_MIN_THRESHOLD, "0");
				//props->add(1, STR_SOLIDITY_KIND, "sum");
				//props->add(1, STR_URI_OUTPUT_TMP, tempDir);
				//cout << tempDir << endl;

				SortingCountAlgorithm<span> algo (filteredBank, *config, repartitor,
						SortingCountAlgorithm<span>::getDefaultProcessorVector (*config, props, solidStorage),
						props);

				algo.execute();

				System::file().rmdir(tempDir);
			}

			//cout << "heo" << endl;
			//delete config;
			//cout << "heo" << endl;
			writeFinishSignal(p);
			//cout << "heo" << endl;
		}

		void writeFinishSignal(Parameter& p){

			string finishFilename = p.outputDir + "/count_synchro/" +  p.bankName + ".ok";
			IFile* file = System::file().newFile(finishFilename, "w");
			delete file;
		}


    };

};

/********************************************************************************/
/*                       Dump solid kmers in ASCII format                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
    	SimkaCount().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
//! [snippet1]
