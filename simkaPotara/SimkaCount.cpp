//! [snippet1]
// We include what we need for the test

#include "SimkaAlgorithm.hpp"
//#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;



template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:

	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads) : BankDelegate (ref), _filter(filter)  {
		//_nbReadsPerDataset = nbReadsPerDataset;
		_maxReads = maxReads;
	}

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

    	//cout << "lala" << endl;
        // We create one iterator from the reference
        Iterator<Sequence>* it = _ref->iterator ();

        // We get the composition for this iterator
        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

		SimkaTruncateIterator<Sequence>* truncIt = new SimkaTruncateIterator<Sequence>(iterators[0], _maxReads+_maxReads/5);
		Filter filter(_filter);
		filter.setMaxReads(_maxReads);
    	FilterIterator<Sequence,Filter>* filterIt = new FilterIterator<Sequence,Filter> (truncIt, filter);
    	return filterIt;

    }

private:

	//vector<u_int64_t> _nbReadsPerDataset;
    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadToProcess;
    size_t _datasetId;
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
        //getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));

        getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);
        if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    }

    void execute ()
    {


    	//size_t datasetId =  getInput()->getInt(STR_ID);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	string bankName =  getInput()->getStr("-bank-name");
    	size_t minReadSize =  getInput()->getInt(STR_SIMKA_MIN_READ_SIZE);
    	double minReadShannonIndex =  getInput()->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
    	u_int64_t maxReads =  getInput()->getInt(STR_SIMKA_MAX_READS);

    	Parameter params(*this, kmerSize, outputDir, bankName, minReadSize, minReadShannonIndex, maxReads);

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
        Parameter (SimkaCount& tool, size_t kmerSize, string outputDir, string bankName, size_t minReadSize, double minReadShannonIndex, u_int64_t maxReads) :
        	tool(tool), kmerSize(kmerSize), outputDir(outputDir), bankName(bankName), minReadSize(minReadSize), minReadShannonIndex(minReadShannonIndex), maxReads(maxReads)  {}
        SimkaCount& tool;
        //size_t datasetId;
        size_t kmerSize;
        string outputDir;
        string bankName;
        size_t minReadSize;
        double minReadShannonIndex;
        u_int64_t maxReads;
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
				IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, p.maxReads);
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
