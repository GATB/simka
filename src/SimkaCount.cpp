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

#include "SimkaPotara.hpp"
#include "minikc/MiniKC.hpp"
//#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

//#define NB_COUNT_CACHE 1
//#define TRACK_DISK_USAGE





template<size_t span=KMER_DEFAULT_SPAN>
class StorageItKmerCount
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t> Kmer_Count;
    //typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    //typedef typename ModelCanonical::Kmer                                   KmerType;

    StorageItKmerCount(Iterator<Count>* it){
    	_it = it;
    	//cout << h5filename << endl;
    	//_bankId = bankId;
    	//_partitionId = partitionId;



		//Iterator<Count>* it2 = partition1.iterator();
		//Collection<Count>& kmers1 = (*partition1)[_partitionId];
		//collections.push_back(&kmers1);

		//_it = kmers1.iterator();

		//_nbKmers = it->estimateNbItems();
		//it2->first();
		//while(!it2->isDone()){
		//	cout << it2->item().value.toString(31) << endl;
		//	it2->next();
		//}
    }

    ~StorageItKmerCount(){
    	delete _it;
    }

    //void setPartitionId(size_t partitionId){
    //	_partitionId = partitionId;
    //}

	bool next(){
		_it->next();

		//cout << "is done?" <<  _it->isDone() << endl;
		return !_it->isDone();
	}

	Count& item(){
		return _it->item();
	}
	//Type& value(){
	//	return _it->item().value;
	//}

	//u_int64_t& abundance(){
	//	return _it->item().abundance;
	//}

	//u_int16_t _bankId;
	//u_int16_t _partitionId;
    Iterator<Count>* _it;
    //u_int64_t _nbKmers;
};



template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:

	Iterator<Sequence>* _it;

	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads, size_t nbDatasets) : BankDelegate (ref), _filter(filter)  {
		//_nbReadsPerDataset = nbReadsPerDataset;
		_maxReads = maxReads;
		_nbDatasets = nbDatasets;
	}


	~SimkaPotaraBankFiltered(){
		delete _it;
	}

    Iterator<Sequence>* iterator ()
    {

        _it = _ref->iterator ();
        //std::vector<Iterator<Sequence>*> iterators = it->getComposition();
        return new SimkaInputIterator<Sequence, Filter> (_it, _nbDatasets, _maxReads, _filter);
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
        getParser()->push_back (new OptionOneParam ("-bank-index",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MAX_READS,   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-datasets",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-partitions",   "bank name", true));
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
    	size_t bankIndex =  getInput()->getInt("-bank-index");
    	size_t minReadSize =  getInput()->getInt(STR_SIMKA_MIN_READ_SIZE);
    	double minReadShannonIndex =  getInput()->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
    	u_int64_t maxReads =  getInput()->getInt(STR_SIMKA_MAX_READS);
    	size_t nbDatasets =   getInput()->getInt("-nb-datasets");
    	size_t nbPartitions =   getInput()->getInt("-nb-partitions");
    	CountNumber abundanceMin =   getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
    	CountNumber abundanceMax =   getInput()->getInt(STR_KMER_ABUNDANCE_MAX);

    	Parameter params(*this, kmerSize, outputDir, bankName, minReadSize, minReadShannonIndex, maxReads, nbDatasets, nbPartitions, abundanceMin, abundanceMax, bankIndex);

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
        Parameter (SimkaCount& tool, size_t kmerSize, string outputDir, string bankName, size_t minReadSize, double minReadShannonIndex, u_int64_t maxReads, size_t nbDatasets, size_t nbPartitions, CountNumber abundanceMin, CountNumber abundanceMax, size_t bankIndex) :
        	tool(tool), kmerSize(kmerSize), outputDir(outputDir), bankName(bankName), minReadSize(minReadSize), minReadShannonIndex(minReadShannonIndex), maxReads(maxReads), nbDatasets(nbDatasets), nbPartitions(nbPartitions), abundanceMin(abundanceMin), abundanceMax(abundanceMax), bankIndex(bankIndex)  {}
        SimkaCount& tool;
        //size_t datasetId;
        size_t kmerSize;
        string outputDir;
        string bankName;
        size_t minReadSize;
        double minReadShannonIndex;
        u_int64_t maxReads;
        size_t nbDatasets;
        size_t nbPartitions;
        CountNumber abundanceMin;
        CountNumber abundanceMax;
        size_t bankIndex;
        size_t localNbPartitions;
    };

    template<size_t span> struct Functor  {

        typedef typename Kmer<span>::Type  Type;
        typedef typename Kmer<span>::Count Count;
        typedef tuple<Count, StorageItKmerCount<span>*> KmerCount_It;
        typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    	struct kxpcomp { bool operator() (KmerCount_It& l,KmerCount_It& r) { return (get<0>(r).value < get<0>(l).value); } } ;

    	void operator ()  (Parameter p){


			IProperties* props = p.tool.getInput();
			vector<string> outInfo;



			IBank* bank = Bank::open(p.outputDir + "/input/" + p.bankName);
			LOCAL(bank);

			/*
			u_int64_t nbSeqs = 1;
	        IBank* sampleBank = new SimkaBankSample(bank, nbSeqs);
			SortingCountAlgorithm<span> sortingCount (sampleBank, props);
			SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();
			sortingCount.addProcessor (proc);
			sortingCount.execute();
			Configuration config = sortingCount.getConfig();
			//_nbPartitions = _maxJobMerge;
			config._nb_partitions = p.nbPartitions;

			uint64_t memoryUsageCachedItems;
			config._nb_cached_items_per_core_per_part = 1 << 8; // cache at least 256 items (128 here, then * 2 in the next while loop)
			do
			{
				config._nb_cached_items_per_core_per_part *= 2;
				memoryUsageCachedItems = 1LL * config._nb_cached_items_per_core_per_part *config._nb_partitions * config._nbCores * sizeof(Type);
			}
			while (memoryUsageCachedItems < config._max_memory * MBYTE / 10);
			*/





			//Configuration config;
			//{
			//Repartitor* repartitor = new Repartitor();
			//LOCAL(repartitor);

			//{
			//	Storage* storage = StorageFactory(STORAGE_HDF5).load (p.outputDir + "/" + "config.h5");
					//	LOCAL (storage);
			//	config.load(storage->getGroup(""));
			//	repartitor->load(storage->getGroup(""));
			//}

				//config._abundanceUserNb = 1;
				//config._abundance.clear();
				//CountRange range(props->getInt(STR_KMER_ABUNDANCE_MIN), 100000);
				//config._abundance.push_back(range);

				/*
				vector<size_t> cacheIndexes;
				cacheIndexes.resize(p.nbPartitions);
				vector<vector<Count> > caches;
	        	caches.resize(p.nbPartitions);
		    	for(size_t i=0; i<p.nbPartitions; i++){
		    		caches[i].resize(NB_COUNT_CACHE);
		    		cacheIndexes[i] = 0;
		    	}
				 */

				//string outputDir = p.outputDir + "/solid/" + p.bankName;
				//System::file().mkdir(outputDir, -1);

			u_int64_t nbReads = 0;
			string tempDir = p.outputDir + "/temp/" + p.bankName;
			System::file().mkdir(tempDir, -1);

			{

				//cout << i << endl;
				//string outputDir = p.outputDir + "/comp_part" + to_string(p.datasetId) + "/";

				//cout << "\tinput: " << p.outputDir + "/input/" + p.bankName << endl;

				SimkaSequenceFilter sequenceFilter(p.minReadSize, p.minReadShannonIndex);
				IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, p.maxReads, p.nbDatasets);
				// = new SimkaPotaraBankFiltered(bank)
				LOCAL(filteredBank);
				//LOCAL(bank);


				string solidsName = p.outputDir + "/solid/" +  p.bankName + ".h5";
				Storage* solidStorage = StorageFactory(STORAGE_HDF5).create (solidsName, true, false);

	    		ConfigurationAlgorithm<span> configAlgo(filteredBank, props);
	    		configAlgo.execute();
	            RepartitorAlgorithm<span> repart (filteredBank, solidStorage->getGroup(""), configAlgo.getConfiguration());
	            repart.execute ();
	            Repartitor* repartitor = new Repartitor();
	            LOCAL(repartitor);
	            repartitor->load(solidStorage->getGroup(""));


				SimkaCompressedProcessor<span>* proc = new SimkaCompressedProcessor<span>(p.abundanceMin, p.abundanceMax, p.bankIndex);
				CountProcessorDump<span>* procDump = new CountProcessorDump     <span> (solidStorage->getGroup("dsk"), p.kmerSize);




				if(p.kmerSize <= 15){
					cout << "Mini Kc a remettre" << endl;
					//MiniKC<span> miniKc(p.tool.getInput(), p.kmerSize, filteredBank, *repartitor, proc);
					//miniKc.execute();

					//nbReads = miniKc._nbReads;
				}
				else{
					//SimkaCompressedProcessor<span>* proc = new SimkaCompressedProcessor<span>(bags, caches, cacheIndexes, p.abundanceMin, p.abundanceMax);
					std::vector<ICountProcessor<span>* > procs;
					//procs.push_back(proc);
					ICountProcessor<span>* result = 0;

					result = new CountProcessorChain<span> (
							proc,
							procDump,
							NULL
					);

					/** We set some name. */
					result->setName ("dsk");

					//SortingCountAlgorithm<span> algo (filteredBank, props);
					SortingCountAlgorithm<span> algo (filteredBank, configAlgo.getConfiguration(), repartitor,
							procs,
							props);
					algo.addProcessor(result);
					algo.execute();

					//delete proc;
					//delete procDump;
					//delete result;
					nbReads = algo.getInfo()->getInt("seq_number");
					p.localNbPartitions = algo.getConfig()._nb_partitions;
				}
			}

			dispatchKmerCounts(p);

			u_int64_t nbDistinctKmers = 0;
			u_int64_t nbKmers = 0;
			u_int64_t chord_N2 = 0;
			for(size_t i=0; i<p.nbPartitions; i++){
				nbDistinctKmers += _nbDistinctKmerPerParts[i];
				nbKmers += _nbKmerPerParts[i];
				chord_N2 += _chordNiPerParts[i];
			}
			//cout << nbDistinctKmers << endl;

			//cout << "CHECK NB READS PER DATASET:  " << nbReads << endl;
			outInfo.push_back(Stringify::format("%llu", nbReads));
			outInfo.push_back(Stringify::format("%llu", nbDistinctKmers));
			outInfo.push_back(Stringify::format("%llu", nbKmers));
			outInfo.push_back(Stringify::format("%llu", chord_N2));

			string contents = "";
			for(size_t i=0; i<_nbDistinctKmerPerParts.size(); i++){
				contents += Stringify::format("%llu", _nbDistinctKmerPerParts[i]) + "\n";
			}
			IFile* nbKmerPerPartFile = System::file().newFile(p.outputDir + "/kmercount_per_partition/" + p.bankName + ".txt", "w");
			nbKmerPerPartFile->fwrite(contents.c_str(), contents.size(), 1);
			nbKmerPerPartFile->flush();
			delete nbKmerPerPartFile;

#ifdef TRACK_DISK_USAGE
				string command = "du -sh " +  p.outputDir;
				system(command.c_str());
#endif

				System::file().rmdir(tempDir);



		    	//delete proc;
				//	}



			//cout << "heo" << endl;
			//delete config;
			//cout << "heo" << endl;
			writeFinishSignal(p, outInfo);
			//cout << "heo" << endl;
		}

    	vector<u_int64_t> _nbKmerPerParts;
    	vector<u_int64_t> _nbDistinctKmerPerParts;
    	vector<u_int64_t> _chordNiPerParts;

    	void dispatchKmerCounts(Parameter p){

    		_nbKmerPerParts = vector<u_int64_t>(p.nbPartitions, 0);
    		_nbDistinctKmerPerParts =  vector<u_int64_t>(p.nbPartitions, 0);
    		_chordNiPerParts = vector<u_int64_t>(p.nbPartitions, 0);


			vector<Bag<Kmer_BankId_Count>* > bags;
			vector<Bag<Kmer_BankId_Count>* > cachedBags;
	    	for(size_t i=0; i<p.nbPartitions; i++){
				string outputFilename = p.outputDir + "/solid/part_" + Stringify::format("%i", i) + "/__p__" + Stringify::format("%i", p.bankIndex) + ".gz";
				Bag<Kmer_BankId_Count>* bag = new BagGzFile<Kmer_BankId_Count>(outputFilename);
				Bag<Kmer_BankId_Count>* cachedBag = new BagCache<Kmer_BankId_Count>(bag, 10000);
				cachedBags.push_back(cachedBag);
				//BagCache bagCache(*bag, 10000);
	        	bags.push_back(bag);
	    	}


			Storage* solidStorage = 0;
			string solidsName = p.outputDir + "/solid/" +  p.bankName + ".h5";
			//cout << solidsName << endl;
			//bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
			solidStorage = StorageFactory(STORAGE_HDF5).load(solidsName);
			LOCAL(solidStorage);

			vector<StorageItKmerCount<span>*> its;
			Partition<Count>& solidKmers = solidStorage->getGroup("dsk").getPartition<Count>("solid");
			for(size_t i=0; i<p.localNbPartitions; i++){
				//cout << "lala" << endl;
				//cout << solidKmers[i].iterable()->getNbItems () << endl;
				its.push_back(new StorageItKmerCount<span>(solidKmers[i].iterable()->iterator()));
			}

			std::priority_queue< KmerCount_It, vector<KmerCount_It>,kxpcomp > pq;
			//StorageItKmerCount<span>* bestIt;
			size_t bestPart;

			for(size_t i=0; i<its.size(); i++){
				its[i]->_it->first();
				//StorageItKmerCount<span>* it = its[i];
				//it->_it->first();
			}

			//fill the  priority queue with the first elems
			for (size_t ii=0; ii<its.size(); ii++)
			{
				//pq.push(Kmer_BankId_Count(ii,its[ii]->value()));
				pq.push(KmerCount_It(its[ii]->item(), its[ii]));
			}

			StorageItKmerCount<span>* bestIt;

			if (pq.size() != 0) // everything empty, no kmer at all
			{
				//get first pointer
				//bestPart =
				//bestIt = get<3>(pq.top()); pq.pop();
				KmerCount_It kmerCountIt = pq.top(); pq.pop();
				size_t part = oahash(get<0>(kmerCountIt).value) % p.nbPartitions;
				cachedBags[part]->insert(Kmer_BankId_Count(get<0>(kmerCountIt).value, p.bankIndex, get<0>(kmerCountIt).abundance));
				_nbDistinctKmerPerParts[part] += 1;
				_nbKmerPerParts[part] += get<0>(kmerCountIt).abundance;
				_chordNiPerParts[part] += pow(get<0>(kmerCountIt).abundance, 2);

				bestIt = get<1>(kmerCountIt);


				while(1){

					if (! bestIt->next())
					{
						//reaches end of one array
						if(pq.size() == 0){
							break;
						}

						//otherwise get new best
						//best_p = get<1>(pq.top()) ; pq.pop();
						bestIt = get<1>(pq.top()); pq.pop();
					}

					pq.push(KmerCount_It(bestIt->item(), bestIt));
					//pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

			    	bestIt = get<1>(pq.top()); pq.pop();
			    	//_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
					size_t part = oahash(bestIt->item().value) % p.nbPartitions;
					//cout << part << " " << p.nbPartitions << endl;
					//cout << bestIt->item().value.toString(21) << " " << bestIt->item().abundance << endl;
					cachedBags[part]->insert(Kmer_BankId_Count(bestIt->item().value, p.bankIndex, bestIt->item().abundance));

					//Kmer_BankId_Count item(kmer, _bankIndex, count[0]);
					_nbDistinctKmerPerParts[part] += 1;
					_nbKmerPerParts[part] += bestIt->item().abundance;
					_chordNiPerParts[part] += pow(bestIt->item().abundance, 2);
				}
			}


	    	for(size_t i=0; i<p.nbPartitions; i++){
	    		//bags[i]->flush();
	    		//cachedBags[i]->flush();
	    		delete cachedBags[i];
	    		//delete bags[i];
	    	}

    	}


		void writeFinishSignal(Parameter& p, const vector<string>& outInfo){

			string finishFilename = p.outputDir + "/count_synchro/" +  p.bankName + ".ok";
			IFile* file = System::file().newFile(finishFilename, "w");
			string contents = "";

			for(size_t i=0; i<outInfo.size(); i++){
				contents += outInfo[i] + "\n";
			}
			file->fwrite(contents.c_str(), contents.size(), 1);
			file->flush();

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
