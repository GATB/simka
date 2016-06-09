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


#include <gatb/gatb_core.hpp>
#include <SimkaAlgorithm.hpp>
#include <SimkaDistance.hpp>

#include "../thirdparty/IteratorKmerH5/IteratorKmerH5.hpp"
#include "../thirdparty/quasi_dictionary/src/quasidictionary.h"

// We use the required packages
using namespace std;







struct Parameter
{
    Parameter (IProperties* props, string inputFilename, string outputDir, size_t partitionId, size_t kmerSize, double minShannonIndex, bool computeSimpleDistances, bool computeComplexDistances, size_t nbPartitions, string inputQueryFilename) : props(props), inputFilename(inputFilename), outputDir(outputDir), partitionId(partitionId), kmerSize(kmerSize), minShannonIndex(minShannonIndex), computeSimpleDistances(computeSimpleDistances), computeComplexDistances(computeComplexDistances), nbPartitions(nbPartitions), inputQueryFilename(inputQueryFilename) {}
    IProperties* props;
    string inputFilename;
    string outputDir;
    size_t partitionId;
    size_t kmerSize;
    double minShannonIndex;
    bool computeSimpleDistances;
    bool computeComplexDistances;
    size_t nbPartitions;
    string inputQueryFilename;
};


template<size_t span=KMER_DEFAULT_SPAN>
class StorageIt
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    //typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    //typedef typename ModelCanonical::Kmer                                   KmerType;

    StorageIt(Iterator<Count>* it, size_t bankId, size_t partitionId){
    	_it = it;
    	//cout << h5filename << endl;
    	_bankId = bankId;
    	_partitionId = partitionId;



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

    ~StorageIt(){

    }

    //void setPartitionId(size_t partitionId){
    //	_partitionId = partitionId;
    //}

	bool next(){
		_it->next();

		//cout << "is done?" <<  _it->isDone() << endl;
		return !_it->isDone();
	}

	Type& value(){
		return _it->item().value;
	}

	CountNumber& abundance(){
		return _it->item().abundance;
	}

	u_int16_t getBankId(){
		return _bankId;
	}

	//u_int64_t getNbKmers(){
	//	return _nbKmers;
	//}

	u_int16_t _bankId;
	u_int16_t _partitionId;
    Iterator<Count>* _it;
    //u_int64_t _nbKmers;
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
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank, CountNumber abundance)  {  _abundancePerBank [idxBank] += abundance;  }

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

private:
    CountVector& _abundancePerBank;
};





template<size_t span>
class SimkaMergeAlgorithm : public Algorithm
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;

	typedef std::pair<u_int16_t, Type> kxp; //id pointer in vec_pointer , value
	struct kxpcomp { bool operator() (kxp l,kxp r) { return ((r.second) < (l.second)); } } ;

	Parameter& p;

	SimkaMergeAlgorithm(Parameter& p) :
		Algorithm("SimkaMergeAlgorithm", 1, p.props), p(p)
	{
		_abundanceThreshold.first = 0;
		_abundanceThreshold.second = 999999999;
	}

	~SimkaMergeAlgorithm(){
		delete _progress;
	}

	//pthread_t statThread;_datasetNbReads

	void createInfo(Parameter& p){

    	for(size_t i=0; i<_nbBanks; i++){
    		string name = _datasetIds[i];
    		string countFilename = p.outputDir + "/count_synchro/" +  name + ".ok";

    		string line;
	    	ifstream file(countFilename.c_str());
	    	vector<string> lines;
			while(getline(file, line)){
				if(line == "") continue;
				lines.push_back(line);
			}
			file.close();

			u_int64_t nbReads = strtoull(lines[0].c_str(), NULL, 10);
			_stats->_datasetNbReads.push_back(nbReads);
			_stats->_nbSolidDistinctKmersPerBank[i] = strtoull(lines[1].c_str(), NULL, 10);
			_stats->_nbSolidKmersPerBank[i] = strtoull(lines[2].c_str(), NULL, 10);
			_stats->_chord_sqrt_N2[i] = sqrt(strtoull(lines[3].c_str(), NULL, 10));
    	}

	}

	void execute(){

		_isQuery = p.inputQueryFilename != "__0__";
		removeStorage(p);

		_partitionId = p.partitionId;

		createDatasetIdList(p);
		_nbBanks = _datasetIds.size();

		//SimkaDistanceParam distanceParams(p.props);
		_stats = new SimkaStatistics(_nbBanks, p.computeSimpleDistances, p.computeComplexDistances);

		createInfo(p);

		if(_isQuery){
			createQuery(p);
		}
		//createBloom(p);

		//createProcessor(p);
		_processor = new SimkaCountProcessorSimple<span> (_stats, _nbBanks, p.kmerSize, _abundanceThreshold, SUM, false, p.minShannonIndex);
		//_processor->use();

		vector<StorageIt<span>*> its;
		//vector<Partition<Count>*> partitions;
		//vector<Collection<Count>*> collections;
		//vector<Iterator<Count>*> its;
		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;
		//vector<Storage*> storages;

		//size_t nbPartitions;
		u_int64_t nbKmers = 0;
		u_int64_t nbKmersProcessed = 0;
		string line;


		vector<IterableGzFile<Count>* > partitions;
		//vector<Iterator<Count>* > partitionIts;
    	for(size_t i=0; i<_nbBanks; i++){
    		string filename = p.outputDir + "/solid/" +  _datasetIds[i] + "/" + "part" + Stringify::format("%i", _partitionId);
    		//cout << filename << endl;
    		IterableGzFile<Count>* partition = new IterableGzFile<Count>(filename);
    		partitions.push_back(partition);
    		its.push_back(new StorageIt<span>(partition->iterator(), i, _partitionId));
    		//nbKmers += partition->estimateNbItems();

    		size_t currentPart = 0;
	    	ifstream file((p.outputDir + "/kmercount_per_partition/" +  _datasetIds[i] + ".txt").c_str());
			while(getline(file, line)){
				if(line == "") continue;
				if(currentPart == _partitionId){
					//cout << stoull(line) << endl;
					nbKmers += strtoull(line.c_str(), NULL, 10);
					break;
				}
				currentPart += 1;
			}
			file.close();
    	}



		//vector<Partition<Count>* > partitions;
		//for(size_t i=0; i<_nbBanks; i++){

			/*
			string solidH5Filename = p.outputDir + "/solid/" +  _datasetIds[i] + ".h5";
			Storage* storage1 = StorageFactory(STORAGE_HDF5).load (solidH5Filename);
			Group& dskGroup1 = storage1->root().getGroup("dsk");
			string nbPartitionsStrg = dskGroup1.getGroup("solid").getProperty ("nb_partitions");
			nbPartitions = atol (nbPartitionsStrg.c_str());
			Partition<Count>* partition1 = &dskGroup1.getPartition<Count>("solid");

			partitions.push_back(partition1);*/
		//}

    	/*
		for(size_t i=0; i<_nbBanks; i++){

			//string solidH5Filename = p.outputDir + "/solid/" +  _datasetIds[i] + ".h5";
			StorageIt<span>* it = new StorageIt<span>(partitions[i], i, p.partitionId, nbPartitions);
			its.push_back(it);
			nbKmers += it->getNbKmers();

			//partitions[i]->remove();
			//partitions[i]->forget();
		}*/

		//partitions.clear();

    	u_int64_t progressStep = nbKmers / 1000;
		_progress = new ProgressSynchro (
			createIteratorListener (nbKmers, "Merging kmers"),
			System::thread().newSynchronizer());
		_progress->init ();







		for(size_t i=0; i<_nbBanks; i++){

			StorageIt<span>* it = its[i];
			it->_it->first();

			//partitionIts[i]->first();

			//while(!it->_it->isDone()){

			//	it->_it->next();
			//	cout << it->_it->item().value.toString(_kmerSize) << " " << it->_it->item().abundance << endl;
			//}

		}


		u_int16_t best_p;
		Type previous_kmer;

	    CountVector abundancePerBank(_nbBanks, 0);
		SimkaCounterBuilderMerge solidCounter(abundancePerBank);
		size_t nbBankThatHaveKmer = 0;

	    //fill the  priority queue with the first elems
	    for (size_t ii=0; ii<_nbBanks; ii++)
	    {
	        //if(its[ii]->next())  {  pq.push(kxp(ii,its[ii]->value()));  }
	    	pq.push(kxp(ii,its[ii]->value()));
	    }

    	//cout << "lala " << pq.size() << endl;

	    if (pq.size() != 0) // everything empty, no kmer at all
	    {
	        //get first pointer
	        best_p = pq.top().first ; pq.pop();

	        previous_kmer = its[best_p]->value();

	        solidCounter.init (its[best_p]->getBankId(), its[best_p]->abundance());
	        nbBankThatHaveKmer = 1;

	        //merge-scan all 'virtual' arrays and output counts
	        while (1)
	        {
	            //go forward in this array or in new array of reaches end of this one
	            if (! its[best_p]->next())
	            {
	                //reaches end of one array
	                if(pq.size() == 0) break; //everything done

	                //otherwise get new best
	                best_p = pq.top().first ; pq.pop();
	            }

	            if (its[best_p]->value() != previous_kmer )
	            {
	                //if diff, changes to new array, get new min pointer
	                pq.push(kxp(best_p,its[best_p]->value())); //push new val of this pointer in pq, will be counted later

	                best_p = pq.top().first ; pq.pop();

	                //if new best is diff, this is the end of this kmer
	                if(its[best_p]->value()!=previous_kmer )
	                {

						nbKmersProcessed += nbBankThatHaveKmer;
						if(nbKmersProcessed > progressStep){
							//cout << "queue size:   " << pq.size() << endl;
							//cout << nbKmersProcessed << endl;
							_progress->inc(nbKmersProcessed);
							nbKmersProcessed = 0;
						}

	                    insert (previous_kmer, abundancePerBank, nbBankThatHaveKmer);

	                    solidCounter.init (its[best_p]->getBankId(), its[best_p]->abundance());
	                    nbBankThatHaveKmer = 1;
	                    previous_kmer = its[best_p]->value();
	                }
	                else
	                {
	                    solidCounter.increase (its[best_p]->getBankId(), its[best_p]->abundance());
	                    nbBankThatHaveKmer += 1;
	                }
	            }
	            else
	            {
	                solidCounter.increase (its[best_p]->getBankId(), its[best_p]->abundance());
	                nbBankThatHaveKmer += 1;
	            }
	        }

        	//cout << nbBankThatHaveKmer << endl;

	    	//cout << previous_kmer.toString(p.kmerSize) << endl;
	        //for(size_t i=0; i<abundancePerBank.size(); i++){
	        //	cout << abundancePerBank[i] << " ";
	        //}
	        //cout << endl;

	        //last elem
            insert (previous_kmer, abundancePerBank, nbBankThatHaveKmer);
	        //if(nbBankThatHaveKmer > 1){
	        //	_processor->process (_partitionId, previous_kmer, abundancePerBank);
	        //}
	        //this->insert (previous_kmer, solidCounter);
	    }

		string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";
		_stats->save(filename); //storage->getGroup(""));


		//saveStats(p);
		writeFinishSignal(p);
		_progress->finish();


		/*
		cout << endl << endl;
		cout << "Reference Breadth: " << endl;
		for(size_t i=0; i<_queriesPresenceAbsences.size(); i++){
			cout << "\t" + _datasetIds[i] + ":  " << (float) _queriesPresenceAbsences[i] / 4564814 << endl;
		}

		cout << endl << "Kmers's abundance" << endl;
		for(size_t i=0; i<_queriesAbundances.size(); i++){
			cout << "\t" + _datasetIds[i] + ":  " << _queriesAbundances[i]<< endl;
		}

		cout << endl;
		cout << "Reference depth" << endl;
		for(size_t i=0; i<_queriesAbundances.size(); i++){
			cout << "\t" + _datasetIds[i] + ":  " << (float) _queriesAbundances[i] / (float) _queriesPresenceAbsences[i] << endl;
		}
		cout << endl;*/
	}


	void createDatasetIdList(Parameter& p){

		string datasetIdFilename = p.outputDir + "/" + "datasetIds";
		IFile* inputFile = System::file().newFile(datasetIdFilename, "rb");
		//IFile* bankFile = System::file().newFile(_banksInputFilename, "wb");

		inputFile->seeko(0, SEEK_END);
		u_int64_t size = inputFile->tell();
		inputFile->seeko(0, SEEK_SET);
		char buffer2[size];
		inputFile->fread(buffer2, size, size);
		string fileContents(buffer2, size);

		string line;
		string linePart;
		vector<string> linePartList;
		stringstream fileContentsStream(fileContents);

		//string bankFileContents = "";

		//u_int64_t lineIndex = 0;

		while(getline(fileContentsStream, line)){

			if(line == "") continue;

			_datasetIds.push_back(line);
		}

		//bankFileContents.erase(bankFileContents.size()-1);
		//bankFileContents.pop_back(); // "remove last /n

		//bankFile->fwrite(bankFileContents.c_str(), bankFileContents.size(), 1);

		delete inputFile;
	}

	void createProcessor(Parameter& p){



		//ICountProcessor<span>* proc = _processor->clone();
		//proc->use();

		//_processors.push_back(proc);
	}

	void insert(const Type& kmer, const CountVector& abundancePerBank, size_t nbBankThatHaveKmer){

		if(_isQuery){
			bool exists;
			_queryKmers->get_value(kmer.getVal(), exists, _queryIds);

			if(exists){

				for(size_t i=0; i<abundancePerBank.size(); i++){
					CountNumber abundance = abundancePerBank[i];
					if(abundance == 0) continue;

					for(size_t j=0; j<_queryIds.size(); j++){
						u_int32_t queryId = _queryIds[j];
						_queriesAbundances[queryId][i] += abundance;
						_queriesPresenceAbsences[queryId][i] += 1;
					}
				}
			}
			//cout << exists << endl;
		}


		if(nbBankThatHaveKmer > 1)
			_processor->process (_partitionId, kmer, abundancePerBank);


	}

	void removeStorage(Parameter& p){
		//Storage* storage = 0;
		//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, true);
		//LOCAL (storage);
	}

	void saveStats(Parameter& p){

		string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";
		//_processor->finishClones(_processors);
		//Storage* storage = 0;
		//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, false);
		//LOCAL (storage);
		_stats->save(filename); //storage->getGroup(""));

		//cout << _stats->_nbKmers << endl;

		//_processors[0]->forget();
		//_processor->forget();

	}

	void writeFinishSignal(Parameter& p){

		if(_isQuery){
			writeQueryResults(p);
		}

		string finishFilename = p.outputDir + "/merge_synchro/" +  SimkaAlgorithm<>::toString(p.partitionId) + ".ok";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
	}

	void writeQueryResults(Parameter& p){

		string filename = p.outputDir + "/merge_synchro/" +  SimkaAlgorithm<>::toString(p.partitionId) + ".query";

		BagGzFile<u_int32_t>* file = new BagGzFile<u_int32_t>(filename);

	    for(size_t i=0; i<_nbQueries; i++){
			for(size_t j=0; j<_nbBanks; j++){
				file->insert(_queriesAbundances[i][j]);
			}
	    }

	    for(size_t i=0; i<_nbQueries; i++){
			for(size_t j=0; j<_nbBanks; j++){
				file->insert(_queriesPresenceAbsences[i][j]);
			}
	    }

		file->flush();

		delete file;
	}

	/*
	void createBloom(Parameter& p){

		string h5filename = p.outputDir + "/" + "query_kmers.h5";

		Storage* storage = StorageFactory(STORAGE_HDF5).load (h5filename);
		LOCAL (storage);

		Partition<Count>& solidCollection = storage->root().getGroup("dsk").getPartition<Count> ("solid");

		Collection<Count>& collection  = solidCollection[p.partitionId];

		//u_int64_t solidFileSize = collection.iterable().getNbItems();
		u_int64_t nb_kmers_infile = collection.iterable()->getNbItems();


		size_t NBITS_PER_KMER = 12;
		u_int64_t estimatedBloomSize = (u_int64_t) ((double)nb_kmers_infile * NBITS_PER_KMER);
		if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

		Iterator<Count>* itKmers = createIterator<Count> (
				collection.iterable()->iterator(),
																		nb_kmers_infile,
																		"Indexing queries"
																		);
		LOCAL (itKmers);
		_queryKmerHash = new Hash16<Type, bool>(100);
		bool lala = true;	
		for(itKmers->first(); !itKmers->isDone(); itKmers->next()){
			
			if(lala){
				cout << itKmers->item().value.toString(31) << endl;
				lala = false;
			}
			_queryKmerHash->insert(itKmers->item().value, true);
		}
		cout << _queryKmerHash->size() << endl;
		//BloomBuilder<span> builder (estimatedBloomSize, 7, p.kmerSize, gatb::core::tools::misc::BLOOM_CACHE, 1, 0);
		//_queryKmerBloom = builder.build (itKmers);
	}*/


	void createQuery(Parameter& p){
		indexQueryKmers(p);

		_queriesAbundances.resize(_nbQueries);
		for(size_t i=0; i<_nbQueries; i++){
			_queriesAbundances[i] = vector<u_int32_t>(_nbBanks, 0);
		}


		_queriesPresenceAbsences.resize(_nbQueries);
		for(size_t i=0; i<_nbQueries; i++){
			_queriesPresenceAbsences[i] = vector<u_int32_t>(_nbBanks, 0);
		}

	}

	void indexQueryKmers(Parameter& p){

		string h5filename = p.outputDir + "/" + "query_kmers.h5";

		Storage* storage = StorageFactory(STORAGE_HDF5).load (h5filename);
		LOCAL (storage);

		Partition<Count>& solidCollection = storage->root().getGroup("dsk").getPartition<Count> ("solid");

		Collection<Count>& collection  = solidCollection[p.partitionId];

		//u_int64_t solidFileSize = collection.iterable().getNbItems();
		u_int64_t nbSolidKmers = collection.iterable()->getNbItems();


		//size_t NBITS_PER_KMER = 12;
		//u_int64_t estimatedBloomSize = (u_int64_t) ((double)nb_kmers_infile * NBITS_PER_KMER);
		//if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

		//Iterator<Count>* itKmers = createIterator<Count> (
		//		collection.iterable()->iterator(),
		//		nbSolidKmers,
		//		"Indexing queries"
		//);
		Iterator<Count>* itKmers = collection.iterable()->iterator();
		LOCAL (itKmers);


		if(nbSolidKmers==0){
			cout<<"No solid kmers in bank -- exit"<<endl;
			exit(0);
		}
		IteratorKmerH5Wrapper<span> iteratorOnKmers (itKmers);
		_queryKmers = new quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper<span>, u_int32_t> (nbSolidKmers, iteratorOnKmers, 8, 2);

		fillQueryKmers(p);
	}

	void fillQueryKmers(Parameter& p){

		IBank* queryBank = Bank::open(p.inputQueryFilename);
		LOCAL(queryBank);

		//Iterator<Sequence>* itSeq = queryBank->iterator();
		Iterator<Sequence>* itSeq = createIterator(queryBank->iterator(), queryBank->estimateNbItems(), "Filling index");
		LOCAL(itSeq);

		//Model definition of a kmer iterator (this one put kmer in cannonical form)
		Kmer<>::ModelCanonical model(p.kmerSize);
		Kmer<>::ModelCanonical::Iterator kmerIt(model);

		Sequence* sequence;
		_nbQueries = 0;

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){
			sequence = &itSeq->item();
			kmerIt.setData (sequence->getData());

			for (kmerIt.first(); !kmerIt.isDone(); kmerIt.next()){

				//cout << _kmerIt->value().toString(kmerSize) << endl;

				u_int64_t kmer = kmerIt->value().getVal();

				u_int32_t value = sequence->getIndex();
				_queryKmers->set_value(kmer, value);
			}

			_nbQueries += 1;
		}
	}







private:
	size_t _nbBanks;
	pair<size_t, size_t> _abundanceThreshold;
	vector<string> _datasetIds;
	size_t _partitionId;
	SimkaStatistics* _stats;
	SimkaCountProcessorSimple<span>* _processor;
	//vector<ICountProcessor<span>*> _processors;

	IteratorListener* _progress;
	IBloom<Type>* _queryKmerBloom;
	Hash16<Type, bool>* _queryKmerHash;
	quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper<span>, u_int32_t>* _queryKmers;

	vector<vector<u_int32_t>> _queriesAbundances;
	vector<vector<u_int32_t>> _queriesPresenceAbsences;
	bool _isQuery;
	u_int64_t _nbQueries;
	vector<u_int32_t> _queryIds;
	//vector<u_int64_t> queryAbundances(_nbBanks, 0);
	//vector<u_int64_t> queryAbundancesCoverages(_nbBanks, 0);
};


class SimkaMerge : public Tool
{
public:

	SimkaMerge () : Tool ("SimkaMerge")
    {
		//Original input filename given to simka. Used to recreate dataset id list
        getParser()->push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size", true));
        getParser()->push_back (new OptionOneParam (STR_URI_INPUT,   "input filename", true));
        getParser()->push_back (new OptionOneParam ("-out-tmp-simka",   "tmp output", true));
        getParser()->push_back (new OptionOneParam ("-partition-id",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX,   "bank name", true));

        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES.c_str(), "compute simple distances"));
        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES.c_str(), "compute complex distances"));
        getParser()->push_back (new OptionOneParam ("-nb-partitions",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-in-query",   "bank name", true));
    }

    void execute ()
    {


    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	size_t partitionId =  getInput()->getInt("-partition-id");
    	string inputFilename =  getInput()->getStr(STR_URI_INPUT);
    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	double minShannonIndex =   getInput()->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
    	bool computeSimpleDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
    	bool computeComplexDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);
    	size_t nbPartitions =   getInput()->getInt("-nb-partitions");
    	string inputQueryFilename =   getInput()->getStr("-in-query");

    	Parameter params(getInput(), inputFilename, outputDir, partitionId, kmerSize, minShannonIndex, computeSimpleDistances, computeComplexDistances, nbPartitions, inputQueryFilename);

        Integer::apply<Functor,Parameter> (kmerSize, params);

    }




    template<size_t span>
    struct Functor  {

    	void operator ()  (Parameter& p)
		{
    		SimkaMergeAlgorithm<span>(p).execute();
		}

    };
};


int main (int argc, char* argv[])
{
    try
    {
    	SimkaMerge().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}


//! [snippet1]
