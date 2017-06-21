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

// We use the required packages
using namespace std;



using namespace gatb::core::system;
using namespace gatb::core::system::impl;


#define MERGE_BUFFER_SIZE 1000
#define SIMKA_MERGE_MAX_FILE_USED 200














struct sortItem_Size_Filename_ID{

	u_int64_t _size;
	size_t _datasetID;

	sortItem_Size_Filename_ID(){}

	sortItem_Size_Filename_ID(u_int64_t size, size_t datasetID){
		_size = size;
		_datasetID = datasetID;
	}
};

bool sortFileBySize (sortItem_Size_Filename_ID i, sortItem_Size_Filename_ID j){
	return ( i._size < j._size );
}

u_int64_t getFileSize(const string& filename){
	std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
	u_int64_t size = in.tellg();
	in.close();
	return size;
}













template<size_t span>
class DistanceCommand : public gatb::core::tools::dp::ICommand //, public gatb::core::system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;


    size_t _bufferIndex;
	size_t _partitionId;
	SimkaStatistics* _stats;
	SimkaCountProcessorSimple<span>* _processor;

	vector<Type> _bufferKmers;
	vector<CountVector> _bufferCounts;

    /** Constructor. */
    DistanceCommand (
    		const string& tmpDir,
			const vector<string>& datasetIds,
    		size_t partitionId,
    		size_t nbBanks,
			bool computeSimpleDistances,
			bool computeComplexDistances,
			size_t kmerSize,
			pair<size_t, size_t>& abundanceThreshold,
			float minShannonIndex
    )
	{
    	_partitionId = partitionId;
		_stats = new SimkaStatistics(nbBanks, computeSimpleDistances, computeComplexDistances, tmpDir, datasetIds);

		_processor = new SimkaCountProcessorSimple<span> (_stats, nbBanks, kmerSize, abundanceThreshold, SUM, false, minShannonIndex);

		_bufferKmers.resize(MERGE_BUFFER_SIZE);
		_bufferCounts.resize(MERGE_BUFFER_SIZE);

		_bufferIndex = 0;
    }

	~DistanceCommand(){
		delete _processor;
		delete _stats;
	}

    //void add(Type& kmer, CountVector& counts){
    //	_bufferIndex +=
    //}
    void setup(size_t bufferIndex, vector<Type>& bufferKmers, vector<CountVector>& bufferCounts){

    	//cout << "hey  " << bufferIndex << endl;
    	_bufferIndex = bufferIndex;

    	for(size_t i=0; i<_bufferIndex; i++){
    		_bufferKmers[i] = bufferKmers[i];
    		_bufferCounts[i] = bufferCounts[i];
    	}
    }

    void execute (){
    	for(size_t i=0; i<_bufferIndex; i++){
    		_processor->process(_partitionId, _bufferKmers[i], _bufferCounts[i]);
    	}
    }

	void use () {}
	void forget () {}
};













struct Parameter
{
    Parameter (IProperties* props, string inputFilename, string outputDir, size_t partitionId, size_t kmerSize, double minShannonIndex, bool computeSimpleDistances, bool computeComplexDistances, size_t nbCores) : props(props), inputFilename(inputFilename), outputDir(outputDir), partitionId(partitionId), kmerSize(kmerSize), minShannonIndex(minShannonIndex), computeSimpleDistances(computeSimpleDistances), computeComplexDistances(computeComplexDistances), nbCores(nbCores) {}
    IProperties* props;
    string inputFilename;
    string outputDir;
    size_t partitionId;
    size_t kmerSize;
    double minShannonIndex;
    bool computeSimpleDistances;
    bool computeComplexDistances;
    size_t nbCores;
};


template<size_t span=KMER_DEFAULT_SPAN>
class StorageIt
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;

    struct Kmer_BankId_Count{
		Type _type;
		u_int32_t _bankId;
		u_int64_t _count;

		Kmer_BankId_Count(){

		}

		Kmer_BankId_Count(Type type, u_int64_t bankId, u_int64_t count){
			_type = type;
			_bankId = bankId;
			_count = count;
		}
	};

    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    //typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    //typedef typename ModelCanonical::Kmer                                   KmerType;

    StorageIt(Iterator<Kmer_BankId_Count>* it, size_t bankId, size_t partitionId){
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

	Type& value(){
		return _it->item()._type;
	}

	u_int16_t getBankId(){
		return _it->item()._bankId;
	}

	u_int64_t& abundance(){
		return _it->item()._count;
	}



	//u_int64_t getNbKmers(){
	//	return _nbKmers;
	//}

	u_int16_t _bankId;
	u_int16_t _partitionId;
    Iterator<Kmer_BankId_Count>* _it;
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













/*
template<size_t span>
class MergeCommand : public gatb::core::tools::dp::ICommand //, public gatb::core::system::SmartPointer
{
public:

        void use () {}
        void forget () {}

    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;


	typedef std::pair<u_int16_t, Type> kxp; //id pointer in vec_pointer , value
	struct kxpcomp { bool operator() (kxp l,kxp r) { return ((r.second) < (l.second)); } } ;

	size_t _currentBuffer;
	u_int64_t _progressStep;
	vector<vector<Type> > _bufferKmers;
	vector<vector<CountVector> > _bufferCounts;
	vector<size_t> _bufferIndex;
	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSharedDistinctKmers;

	MergeCommand (
    		size_t partitionId,
    		size_t nbBanks,
			IteratorListener* progress,
			vector<StorageIt<span>*>& its,
			u_int64_t progressStep,
			size_t nbCores,
			bool computeComplexDistances
    ) :
    	its(its)
	{
		_nbBanks = nbBanks;
    	_partitionId = partitionId;
    	_progress = progress;
    	_progressStep = progressStep;
    	_nbCores = nbCores;
    	_computeComplexDistances = computeComplexDistances;
    	_nbDistinctKmers = 0;
    	_nbSharedDistinctKmers = 0;

		init();
    }

	~MergeCommand(){
		delete solidCounter;	
	}

    //void add(Type& kmer, CountVector& counts){
    //	_bufferIndex +=
    //}
    //void setup(vector<Type>& bufferKmers, vector<CountVector>& bufferCounts){
    //	_bufferKmers = bufferKmers;
    //	_bufferCounts = bufferCounts;
    //}

	size_t _nbCores;
	size_t _partitionId;
	size_t _nbBanks;
	vector<StorageIt<span>*>& its;
	std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;
	u_int64_t nbKmersProcessed;
	IteratorListener* _progress;
	bool _computeComplexDistances;

	u_int16_t best_p;
	Type previous_kmer;
    CountVector abundancePerBank;
	size_t nbBankThatHaveKmer;
	SimkaCounterBuilderMerge* solidCounter;
	bool _isDone;


	void init(){

		_isDone = false;
		solidCounter = new SimkaCounterBuilderMerge(abundancePerBank);
		for(size_t i=0; i<_nbCores; i++){
			vector<Type> vec = vector<Type>(MERGE_BUFFER_SIZE);
			_bufferKmers.push_back(vec);
			vector<CountVector> vec2 = vector<CountVector>(MERGE_BUFFER_SIZE);
			_bufferCounts.push_back(vec2);
			_bufferIndex.push_back(0);
		}


		nbBankThatHaveKmer = 0;
		abundancePerBank.resize(_nbBanks, 0);
		_currentBuffer = 0;
		//_bufferIndex = 0;
		//_bufferSize = 1000;

		nbKmersProcessed = 0;
		//vector<Partition<Count>*> partitions;
		//vector<Collection<Count>*> collections;
		//vector<Iterator<Count>*> its;
		//vector<Storage*> storages;

		//size_t nbPartitions;







		for(size_t i=0; i<_nbBanks; i++){

			StorageIt<span>* it = its[i];
			it->_it->first();

			//partitionIts[i]->first();

			//while(!it->_it->isDone()){

			//	it->_it->next();
			//	cout << it->_it->item().value.toString(_kmerSize) << " " << it->_it->item().abundance << endl;
			//}

		}



	    //fill the  priority queue with the first elems
	    for (size_t ii=0; ii<_nbBanks; ii++)
	    {
	        //if(its[ii]->next())  {  pq.push(kxp(ii,its[ii]->value()));  }
	    	pq.push(kxp(ii,its[ii]->value()));
	    }

	    if (pq.size() != 0) // everything empty, no kmer at all
	    {
	        //get first pointer
	        best_p = pq.top().first ; pq.pop();

	        previous_kmer = its[best_p]->value();

	        solidCounter->init (its[best_p]->getBankId(), its[best_p]->abundance());
	        nbBankThatHaveKmer = 1;
	    }
	}


	void reset(){
		for(size_t i=0; i<_bufferIndex.size(); i++){
			_bufferIndex[i] = 0;
		}
	}


    void execute (){

    	//cout << "lala " << pq.size() << endl;



	        //merge-scan all 'virtual' arrays and output counts
	        while (_currentBuffer < _nbCores)
	        {

	        	//cout << _currentBuffer << endl;
	            //go forward in this array or in new array of reaches end of this one
	            if (! its[best_p]->next())
	            {
	                //reaches end of one array
	                if(pq.size() == 0){
	                	_isDone = true;
	                	break;
	                }

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
						if(nbKmersProcessed > _progressStep){
							//cout << "queue size:   " << pq.size() << endl;
							//cout << nbKmersProcessed << endl;
							_progress->inc(nbKmersProcessed);
							nbKmersProcessed = 0;
						}

						//cout << previous_kmer.toString(p.kmerSize) << endl;
						//for(size_t i=0; i<abundancePerBank.size(); i++){
						//	cout << abundancePerBank[i] << " ";
						//}
						//cout << endl;

				        insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
						//if(nbBankThatHaveKmer > 1)
						//	_processor->process (_partitionId, previous_kmer, abundancePerBank);
	                    //this->insert (previous_kmer, solidCounter);

	                    solidCounter->init (its[best_p]->getBankId(), its[best_p]->abundance());
	                    nbBankThatHaveKmer = 1;
	                    previous_kmer = its[best_p]->value();
	                }
	                else
	                {
	                    solidCounter->increase (its[best_p]->getBankId(), its[best_p]->abundance());
	                    nbBankThatHaveKmer += 1;
	                }
	            }
	            else
	            {
	                solidCounter->increase (its[best_p]->getBankId(), its[best_p]->abundance());
	                nbBankThatHaveKmer += 1;
	            }
	        }

	        if(_isDone){
		        insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
	        }
	        else{
	        }

			_currentBuffer = 0;
			//_bufferIndex = 0;
        	//cout << nbBankThatHaveKmer << endl;

	    	//cout << previous_kmer.toString(p.kmerSize) << endl;
	        //for(size_t i=0; i<abundancePerBank.size(); i++){
	        //	cout << abundancePerBank[i] << " ";
	        //}
	        //cout << endl;

	        //last elem
	        //insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
	        //this->insert (previous_kmer, solidCounter);
	   // }


	    	//cout << "end " << endl;
    }






	void insert(const Type& kmer, const CountVector& counts, size_t nbBankThatHaveKmer){

		_nbDistinctKmers += 1;

        if(_computeComplexDistances || nbBankThatHaveKmer > 1){

        	if(nbBankThatHaveKmer > 1){
        		_nbSharedDistinctKmers += 1;
        	}

			//DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[_currentBuffer]);
			//cmd->_bufferKmers[cmd->_bufferIndex] = kmer;
			//cmd->_bufferCounts[cmd->_bufferIndex] = counts;
        	_bufferKmers[_currentBuffer][_bufferIndex[_currentBuffer]] = kmer;
        	_bufferCounts[_currentBuffer][_bufferIndex[_currentBuffer]] = counts;

			_bufferIndex[_currentBuffer] += 1;
        	if(_bufferIndex[_currentBuffer] >= MERGE_BUFFER_SIZE){
				//DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[_currentBuffer]);
				//cmd->setup(_bufferKmers[_currentBuffer], _bufferCounts[_currentBuffer]);

        		_currentBuffer += 1;
        		if(_currentBuffer >= _nbCores){
        			//dispatch();
        		}
        		else{
        			//_bufferIndex = 0;
        		}
        	}
        	//_processor->process (_partitionId, kmer, counts);
        }

    	//_processor->process (_partitionId, kmer, counts);


		//cout <<_partitiontId << " "<< kmer.toString(31) << endl;
		//_processor->process (_partitionId, kmer, counter.get());
	}


};
*/



template<size_t span>
class DiskBasedMergeSort
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    //typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;

	typedef typename StorageIt<span>::Kmer_BankId_Count Kmer_BankId_Count;

	struct kxp{
		Type _type;
		u_int32_t _bankId;
		u_int64_t _count;
		StorageIt<span>* _it;

		kxp(){

		}

		kxp(Type type, u_int64_t bankId, u_int64_t count, StorageIt<span>* it){
			_type = type;
			_bankId = bankId;
			_count = count;
			_it = it;
		}
	};

	struct kxpcomp { bool operator() (kxp& l, kxp& r) { return (r._type < l._type); } } ;

	string _outputDir;
	string _outputFilename;
	vector<size_t>& _datasetIds;
	size_t _partitionId;
	Bag<Kmer_BankId_Count>* _outputGzFile;
	Bag<Kmer_BankId_Count>* _cachedBag;



    DiskBasedMergeSort(size_t mergeId, const string& outputDir, vector<size_t>& datasetIds, size_t partitionId):
    	_datasetIds(datasetIds)
    {
    	_outputDir = outputDir;
    	_partitionId = partitionId;

    	_outputFilename = _outputDir + "/solid/part_" + Stringify::format("%i", partitionId) + "/__p__" + Stringify::format("%i", mergeId) + ".gz.temp";
    	_outputGzFile = new BagGzFile<Kmer_BankId_Count>(_outputFilename);
    	_cachedBag = new BagCache<Kmer_BankId_Count>(_outputGzFile, 10000);

    }

    ~DiskBasedMergeSort(){
    }

    void execute(){

		vector<IterableGzFile<Kmer_BankId_Count>* > partitions;
		vector<StorageIt<span>*> its;

		size_t _nbBanks = _datasetIds.size();

		for(size_t i=0; i<_nbBanks; i++){
			//cout << _datasetIds[i] << endl;
			string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			//cout << "\t\t" << filename << endl;
			IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);
			partitions.push_back(partition);
			its.push_back(new StorageIt<span>(partition->iterator(), i, _partitionId));
			//nbKmers += partition->estimateNbItems();

			//size_t currentPart = 0;
			//ifstream file((_outputDir + "/kmercount_per_partition/" +  _datasetIds[i] + ".txt").c_str());
			//while(getline(file, line)){
			//	if(line == "") continue;
			//	if(currentPart == _partitionId){
			//		//cout << stoull(line) << endl;
			//		nbKmers += strtoull(line.c_str(), NULL, 10);
			//		break;
			//	}
			//	currentPart += 1;
			//}
			//file.close();
		}

		//u_int64_t progressStep = nbKmers / 1000;
		//_progress = new ProgressSynchro (
		//	createIteratorListener (nbKmers, "Merging kmers"),
		//	System::thread().newSynchronizer());
		//_progress->init ();



		//_nbDistinctKmers = 0;
		//_nbSharedDistinctKmers = 0;
		//u_int64_t nbKmersProcessed = 0;
		//size_t nbBankThatHaveKmer = 0;
		//u_int16_t best_p = 0;
		Type previous_kmer;
		//CountVector abundancePerBank;
		//abundancePerBank.resize(_nbBanks, 0);
		//SimkaCounterBuilderMerge* solidCounter = new SimkaCounterBuilderMerge(abundancePerBank);;
		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;
		StorageIt<span>* bestIt;


		for(size_t i=0; i<_nbBanks; i++){
			StorageIt<span>* it = its[i];
			it->_it->first();
		}

		//fill the  priority queue with the first elems
		for (size_t ii=0; ii<_nbBanks; ii++)
		{
			//pq.push(Kmer_BankId_Count(ii,its[ii]->value()));
			pq.push(kxp(its[ii]->value(), its[ii]->getBankId(), its[ii]->abundance(), its[ii]));
		}

		if (pq.size() != 0) // everything empty, no kmer at all
		{
			//get first pointer
			bestIt = pq.top()._it; pq.pop();
			_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
			//best_p = get<1>(pq.top()) ; pq.pop();
			//previous_kmer = bestIt->value();
			//solidCounter->init (bestIt->getBankId(), bestIt->abundance());
			//nbBankThatHaveKmer = 1;

			while(1){

				if (! bestIt->next())
				{
					//reaches end of one array
					if(pq.size() == 0){
						break;
					}

					//otherwise get new best
					//best_p = get<1>(pq.top()) ; pq.pop();
					bestIt = pq.top()._it; pq.pop();
				}

				pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

		    	bestIt = pq.top()._it; pq.pop();
		    	_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
		    	//cout << bestIt->value().toString(31) << " " << bestIt->getBankId() <<  " "<< bestIt->abundance() << endl;
				//bestIt = get<3>(pq.top()); pq.pop();


				//pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt));

			}


	    	//_outputGzFile->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
	    	//cout << bestIt->value().toString(31) << " " << bestIt->getBankId() <<  " "<< bestIt->abundance() << endl;
		}

		for(size_t i=0; i<partitions.size(); i++){
			delete partitions[i];
		}

		for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}


		_cachedBag->flush();
    	delete _cachedBag;

		for(size_t i=0; i<_nbBanks; i++){
			//cout << _datasetIds[i] << endl;
			string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			System::file().remove(filename);
		}

		string newOutputFilename = _outputFilename;
		newOutputFilename.erase(_outputFilename.size()-5, 5);
    	System::file().rename(_outputFilename, newOutputFilename); //remove .temp at the end of new merged file
    	//_outputFilename = newOutputFilename;
    }

};



template<size_t span>
class SimkaMergeAlgorithm : public Algorithm
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
    //typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;

    //typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;

	typedef typename DiskBasedMergeSort<span>::Kmer_BankId_Count Kmer_BankId_Count;
	typedef typename DiskBasedMergeSort<span>::kxp kxp;


	/*
	struct Kmer_BankId_Count{
		Type _type;
		u_int64_t _bankId;
		u_int64_t _count;

		Kmer_BankId_Count(){

		}

		Kmer_BankId_Count(Type type, u_int64_t bankId, u_int64_t count){
			_type = type;
			_bankId = bankId;
			_count = count;
		}
	};

	struct kxp{
		Type _type;
		u_int32_t _bankId;
		u_int64_t _count;
		StorageIt<span>* _it;

		kxp(){

		}

		kxp(Type type, u_int64_t bankId, u_int64_t count, StorageIt<span>* it){
			_type = type;
			_bankId = bankId;
			_count = count;
			_it = it;
		}
	};*/



	//typedef std::pair<u_int16_t, Type> kxp; //id pointer in vec_pointer , value
    //typedef std::pair<u_int16_t, Type> kxp; //id pointer in vec_pointer , value
	//struct kxpcomp { bool operator() (Kmer_BankId_Count l,Kmer_BankId_Count r) { return ((r.second) < (l.second)); } } ;
	struct kxpcomp { bool operator() (kxp& l,kxp& r) { return (r._type < l._type); } } ;

	Parameter& p;

	SimkaMergeAlgorithm(Parameter& p) :
		Algorithm("SimkaMergeAlgorithm", p.nbCores, p.props), p(p)
	{
		_abundanceThreshold.first = 0;
		_abundanceThreshold.second = 999999999;

		_computeSimpleDistances = p.computeSimpleDistances;
		_computeComplexDistances = p.computeComplexDistances;
		_kmerSize = p.kmerSize;
		_minShannonIndex = p.minShannonIndex;
	}

	~SimkaMergeAlgorithm(){
		//delete _progress;
	}

	//pthread_t statThread;_datasetNbReads

	/*
	void createInfo(Parameter& p){



	}


	void loadCountInfo(){
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

			_stats->_datasetNbReads[i] = nbReads;
			_stats->_nbSolidDistinctKmersPerBank[i] = strtoull(lines[1].c_str(), NULL, 10);
			_stats->_nbSolidKmersPerBank[i] = strtoull(lines[2].c_str(), NULL, 10);
			_stats->_chord_sqrt_N2[i] = sqrt(strtoull(lines[3].c_str(), NULL, 10));
			//cout << _stats->_chord_sqrt_N2[i] << endl;
    	}
	}*/


	//struct sortFileBySize { bool operator() (sortItem_Size_Filename_ID& l,sortItem_Size_Filename_ID& r) { return (r._size < l._size); } } ;

	void execute(){

		_nbCores = p.nbCores;




		removeStorage(p);

		_partitionId = p.partitionId;

		createDatasetIdList(p);
		_nbBanks = _datasetIds.size();

		string partDir = p.outputDir + "/solid/part_" + Stringify::format("%i", _partitionId) + "/";
		vector<string> filenames = System::file().listdir(partDir);
		//cout << filenames.size() << endl;
		vector<string> partFilenames;
		vector<sortItem_Size_Filename_ID> filenameSizes;

		for(size_t i=0; i<filenames.size(); i++){
			if(filenames[i].find("__p__") != std::string::npos){


				string id = string(filenames[i]);
				id.erase(0, 5);
				std::string::size_type pos = id.find(".gz");
				id.erase(pos, 3);

				size_t datasetId = atoll(id.c_str());
				//cout << filenames[i] << " " << datasetId << endl;

				filenameSizes.push_back(sortItem_Size_Filename_ID(getFileSize(partDir+filenames[i]), datasetId));
				//cout << filenames[i] << " " << size << endl;
				//cout << filenames[i] << endl;
			}
		}

		//cout << "mettre un while ici" << endl;
		while(filenameSizes.size() > SIMKA_MERGE_MAX_FILE_USED){

			//cout << "Start merging pass" << endl;
			sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);

			vector<size_t> mergeDatasetIds;
			vector<size_t> toRemoveItem;


			for(size_t i=0; i<SIMKA_MERGE_MAX_FILE_USED; i++){
				sortItem_Size_Filename_ID sfi = filenameSizes[i];
				mergeDatasetIds.push_back(sfi._datasetID);
				//datasetIndex += 1;
				//if(datasetIndex >= _nbBanks) break;

				//cout << mergeDatasetIds[i] << endl;
				//cout << "First val must never be greater than second:   " << i << "  " << _nbBanks << endl;
				//cout << "\t" << get<1>(sfi) << endl;
			}

			for(size_t i=0; i<mergeDatasetIds.size(); i++){
				filenameSizes.erase(filenameSizes.begin());
			}

			size_t mergedId = mergeDatasetIds[0];
			DiskBasedMergeSort<span> diskBasedMergeSort(mergedId, p.outputDir, mergeDatasetIds, _partitionId);
			diskBasedMergeSort.execute();

			filenameSizes.push_back(sortItem_Size_Filename_ID(getFileSize(diskBasedMergeSort._outputFilename), mergedId));

			//cout << "\tmerged id: " <<  mergedId << endl;
			//cout << "\tremainging files: " << filenameSizes.size() << endl;
		}

		//cout << filenameSizes.size() << endl;
		//for(size_t i=0; i<filenameSizes.size(); i++){
		//	cout << filenameSizes[i].first << endl;
		//}

		//size_t nbMerges = 0;
		/*
		//cout << partFilenames.size() << endl;
		exit(1);

		size_t nbMerges = ceil((float)_nbBanks / (float)SIMKA_MERGE_MAX_FILE_USED);
		cout << "nb Merges: " << nbMerges << endl;
		size_t datasetIndex = 0;

		for(size_t i=0; i<nbMerges; i++){

			vector<string> mergeDatasetIds;

			for(size_t j=0; j<SIMKA_MERGE_MAX_FILE_USED; j++){
				mergeDatasetIds.push_back(_datasetIds[datasetIndex]);
				datasetIndex += 1;
				if(datasetIndex >= _nbBanks) break;
			}

			cout << "doivent etre égaux a la dernière passe:    " << _nbBanks << " " << mergeDatasetIds.size() << " " << datasetIndex << endl;

			DiskBasedMergeSort<span> diskBasedMergeSort(i, p.outputDir, mergeDatasetIds, _partitionId);
			diskBasedMergeSort.execute();

		}*/

		//exit(1);
		/* PARALLEL
		for (size_t i=0; i<_nbCores; i++)
	    {
		//cout << i << endl;
	        ICommand* cmd = 0;
	        cmd = new DistanceCommand<span>(p.outputDir, _datasetIds, _partitionId, _nbBanks, _computeSimpleDistances, _computeComplexDistances, _kmerSize, _abundanceThreshold, _minShannonIndex);
	        //cmd->use();
	        _cmds.push_back (cmd);

                        //cout << _cmds[i] << endl;
	    }

		resetCommands();
		*/


		//SimkaDistanceParam distanceParams(p.props);
		//createInfo(p);


		//createProcessor(p);

		//PARALLEL line to remove
		_stats = new SimkaStatistics(_nbBanks, p.computeSimpleDistances, p.computeComplexDistances, p.outputDir, _datasetIds);
		_processor = new SimkaCountProcessorSimple<span> (_stats, _nbBanks, p.kmerSize, _abundanceThreshold, SUM, false, p.minShannonIndex);
		//_processor->use();






		string line;
		vector<IterableGzFile<Kmer_BankId_Count>* > partitions;
		vector<StorageIt<span>*> its;
		u_int64_t nbKmers = 0;

    	for(size_t i=0; i<filenameSizes.size(); i++){
    		size_t datasetId = filenameSizes[i]._datasetID;
    		string filename = p.outputDir + "/solid/part_" + Stringify::format("%i", p.partitionId) + "/__p__" + Stringify::format("%i", datasetId) + ".gz";
    		//cout << filename << endl;
    		IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);
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


		/*
		//vector<Iterator<Count>* > partitionIts;
    	for(size_t i=0; i<_nbBanks; i++){
    		string filename = p.outputDir + "/solid/" +  _datasetIds[i] + "/" + "part" + Stringify::format("%i", _partitionId);
    		//cout << filename << endl;
    		IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 1000);
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
    	}*/

    	//u_int64_t progressStep = nbKmers / 1000;
    	//_progress = new ProgressSynchro (
    	//	createIteratorListener (nbKmers, "Merging kmers"),
    	//	System::thread().newSynchronizer());
    	//_progress->init ();


		/* PARALLEL
		_mergeCommand = new MergeCommand<span>(
    		_partitionId,
    		_nbBanks,
			_progress,
			its,
			progressStep,
			_nbCores,
			p.computeComplexDistances);
		//_mergeCommand->use();
		_cmds.push_back(_mergeCommand);


		//cout << "CMDS SIZE:" << _cmds.size() << endl;


		MergeCommand<span>* mergeCmd = dynamic_cast<MergeCommand<span>*>(_mergeCommand);
		mergeCmd->execute();

		while(!mergeCmd->_isDone){
			//cout << mergeCmd->_isDone << endl;
			//mergeCmd->execute();
			dispatch();
		}

	    dispatch();*/

		_nbDistinctKmers = 0;
		_nbSharedDistinctKmers = 0;
		u_int64_t nbKmersProcessed = 0;
		size_t nbBankThatHaveKmer = 0;
		u_int16_t best_p = 0;
		Type previous_kmer;
	    CountVector abundancePerBank;
		abundancePerBank.resize(_nbBanks, 0);
		SimkaCounterBuilderMerge* solidCounter = new SimkaCounterBuilderMerge(abundancePerBank);;
		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;

    	StorageIt<span>* bestIt;

		for(size_t i=0; i<its.size(); i++){
			StorageIt<span>* it = its[i];
			it->_it->first();
		}

	    //fill the  priority queue with the first elems
	    for (size_t ii=0; ii<its.size(); ii++)
	    {
	    	//pq.push(Kmer_BankId_Count(ii,its[ii]->value()));
	    	pq.push(kxp(its[ii]->value(), its[ii]->getBankId(), its[ii]->abundance(), its[ii]));
	    }

	    if (pq.size() != 0) // everything empty, no kmer at all
	    {
	        //get first pointer
	    	bestIt = pq.top()._it; pq.pop();
	        //best_p = get<1>(pq.top()) ; pq.pop();
	        previous_kmer = bestIt->value();
	        solidCounter->init (bestIt->getBankId(), bestIt->abundance());
	        nbBankThatHaveKmer = 1;

			while(1){

				if (! bestIt->next())
				{
					//reaches end of one array
					if(pq.size() == 0){
						break;
					}

					//otherwise get new best
					//best_p = get<1>(pq.top()) ; pq.pop();
			    	bestIt = pq.top()._it; pq.pop();
				}

		    	//cout << bestIt->value().toString(31) << " " << bestIt->getBankId() <<  " "<< bestIt->abundance() << endl;

				if (bestIt->value() != previous_kmer )
				{
					//if diff, changes to new array, get new min pointer
					pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

			    	bestIt = pq.top()._it; pq.pop();
					//best_p = get<1>(pq.top()) ; pq.pop();

					//if new best is diff, this is the end of this kmer
					if(bestIt->value()!=previous_kmer )
					{

						//nbKmersProcessed += nbBankThatHaveKmer;
						//if(nbKmersProcessed > progressStep){
							//cout << "queue size:   " << pq.size() << endl;
							//cout << nbKmersProcessed << endl;
							//_progress->inc(nbKmersProcessed);
						//nbKmersProcessed = 0;
						//}

						//cout << previous_kmer.toString(p.kmerSize) << endl;
						//for(size_t i=0; i<abundancePerBank.size(); i++){
						//	cout << abundancePerBank[i] << " ";
						//}
						//cout << endl;

						insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
						//if(nbBankThatHaveKmer > 1)
						//	_processor->process (_partitionId, previous_kmer, abundancePerBank);
						//this->insert (previous_kmer, solidCounter);

						solidCounter->init (bestIt->getBankId(), bestIt->abundance());
						nbBankThatHaveKmer = 1;
						previous_kmer = bestIt->value();
					}
					else
					{
						solidCounter->increase (bestIt->getBankId(), bestIt->abundance());
						nbBankThatHaveKmer += 1;
					}
				}
				else
				{
					//cout << "increase" << endl;
					solidCounter->increase (bestIt->getBankId(), bestIt->abundance());
					nbBankThatHaveKmer += 1;
				}
			}

			insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
	    }


		_processor->end();

		//cout << "lala" << endl;
		for(size_t i=0; i<partitions.size(); i++){
			delete partitions[i];
		}


		saveStats(p);


		delete _stats;
		delete _processor;
		//for(size_t i=0; i<its.size(); i++){
		//	delete its[i];
		//}

		/* PARALLEL
		saveStats(p, mergeCmd->_nbDistinctKmers, mergeCmd->_nbSharedDistinctKmers);


		//cout << _cmds.size() << endl;
		for(size_t i=0; i<_cmds.size(); i++){
			//cout << _cmds[i] << endl;
			//_cmds[i]->forget();
			delete _cmds[i];		
		}
		//_cmds.clear();
		//delete _mergeCommand;
		*/
		delete solidCounter;
		for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}

		writeFinishSignal(p);
		//_progress->finish();

	}

	void insert(const Type& kmer, const CountVector& counts, size_t nbBankThatHaveKmer){

		//cout << kmer.toString(31) << endl;
		//for(size_t i=0; i<counts.size(); i++){
		//	cout << counts[i] << " ";
		//}
		//cout << endl;

		_stats->_nbDistinctKmers += 1;

		if(_computeComplexDistances || nbBankThatHaveKmer > 1){

			if(nbBankThatHaveKmer > 1){
				_stats->_nbSharedKmers += 1;
			}

			_processor->process(_partitionId, kmer, counts);

		}
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

	void resetCommands(){
		for (size_t i=0; i<_nbCores; i++){
			DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[i]);
			cmd->_bufferIndex = 0;
		}

	}

	/*
	void dispatch(){


		MergeCommand<span>* mergeCommand = dynamic_cast<MergeCommand<span>*>(_mergeCommand);
		for (size_t i=0; i<_nbCores; i++){
			//cout << mergeCommand->_bufferKmers.size() << endl;
			//cout << i << endl;
			DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[i]);
			cmd->setup(mergeCommand->_bufferIndex[i], mergeCommand->_bufferKmers[i], mergeCommand->_bufferCounts[i]);
		}


		//MergeCommand<span>* mergeCommand = dynamic_cast<MergeCommand<span>*>(_mergeCommand);
		mergeCommand->reset();

		//cout << "start dispatch" << endl;
	    getDispatcher()->dispatchCommands(_cmds, 0);


		//cout << "end dispatch" << endl;
	    resetCommands();


	}*/


	void removeStorage(Parameter& p){
		//Storage* storage = 0;
		//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, true);
		//LOCAL (storage);
	}


	/* PARALLEL
	void saveStats(Parameter& p, const u_int64_t nbDistinctKmers, const u_int64_t nbSharedDistinctKmers){

		_stats = new SimkaStatistics(_nbBanks, p.computeSimpleDistances, p.computeComplexDistances, p.outputDir, _datasetIds);

		for (size_t i=0; i<_nbCores; i++){
			DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[i]);
			cmd->_processor->end();
			(*_stats) += (*cmd->_stats);
		}
		//loadCountInfo();

		string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";

		_stats->_nbDistinctKmers = nbDistinctKmers;
		_stats->_nbSharedKmers = nbSharedDistinctKmers;
		_stats->save(filename); //storage->getGroup(""));

		
		delete _stats;
		//string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";
		//_processor->finishClones(_processors);
		//Storage* storage = 0;
		//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, false);
		//LOCAL (storage);
		//_stats->save(filename); //storage->getGroup(""));

		//cout << _stats->_nbKmers << endl;

		//_processors[0]->forget();
		//_processor->forget();

	}*/

	void saveStats(Parameter& p){

		string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";

		_stats->save(filename); //storage->getGroup(""));


		//string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";
		//_processor->finishClones(_processors);
		//Storage* storage = 0;
		//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, false);
		//LOCAL (storage);
		//_stats->save(filename); //storage->getGroup(""));

		//cout << _stats->_nbKmers << endl;

		//_processors[0]->forget();
		//_processor->forget();

	}

	void writeFinishSignal(Parameter& p){
		string finishFilename = p.outputDir + "/merge_synchro/" +  SimkaAlgorithm<>::toString(p.partitionId) + ".ok";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
	}

private:
	size_t _nbBanks;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	size_t _kmerSize;
	float _minShannonIndex;

	pair<size_t, size_t> _abundanceThreshold;
	vector<string> _datasetIds;
	size_t _partitionId;
	//vector<ICountProcessor<span>*> _processors;

	IteratorListener* _progress;

    vector<ICommand*> _cmds;
	ICommand* _mergeCommand;
	size_t _nbCores;


	SimkaStatistics* _stats;
	SimkaCountProcessorSimple<span>* _processor;
	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSharedDistinctKmers;
};
























class SimkaMerge : public Tool
{
public:

	SimkaMerge () : Tool ("SimkaMerge")
    {
		//Original input filename given to simka. Used to recreate dataset id list
        getParser()->push_back (new OptionOneParam (STR_NB_CORES,   "nb cores", true));
        getParser()->push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size", true));
        getParser()->push_back (new OptionOneParam (STR_URI_INPUT,   "input filename", true));
        getParser()->push_back (new OptionOneParam ("-out-tmp-simka",   "tmp output", true));
        getParser()->push_back (new OptionOneParam ("-partition-id",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX,   "bank name", true));

        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES.c_str(), "compute simple distances"));
        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES.c_str(), "compute complex distances"));
    }

    void execute ()
    {


    	size_t nbCores =  getInput()->getInt(STR_NB_CORES);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	size_t partitionId =  getInput()->getInt("-partition-id");
    	string inputFilename =  getInput()->getStr(STR_URI_INPUT);
    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	double minShannonIndex =   getInput()->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
    	bool computeSimpleDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
    	bool computeComplexDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);

    	Parameter params(getInput(), inputFilename, outputDir, partitionId, kmerSize, minShannonIndex, computeSimpleDistances, computeComplexDistances, nbCores);

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






