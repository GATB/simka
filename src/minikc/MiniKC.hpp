/*
 * MiniKC.hpp
 *
 *  Created on: 16 juin 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_MINIKC_MINIKC_HPP_
#define GATB_SIMKA_SRC_MINIKC_MINIKC_HPP_

#include <gatb/gatb_core.hpp>
//#include "../SimkaCount.cpp"

//typedef u_int16_t CountType;

template<size_t span>
class SimkaCompressedProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;

    //SimkaCompressedProcessor(vector<BagGzFile<Count>* >& bags, vector<vector<Count> >& caches, vector<size_t>& cacheIndexes, CountNumber abundanceMin, CountNumber abundanceMax) : _bags(bags), _caches(caches), _cacheIndexes(cacheIndexes)
    SimkaCompressedProcessor(vector<Bag<Kmer_BankId_Count>* >& bags, vector<u_int64_t>& nbKmerPerParts, vector<u_int64_t>& nbDistinctKmerPerParts, vector<u_int64_t>& chordPerParts, CountNumber abundanceMin, CountNumber abundanceMax, size_t bankIndex) :
    	_bags(bags), _nbDistinctKmerPerParts(nbDistinctKmerPerParts), _nbKmerPerParts(nbKmerPerParts), _chordPerParts(chordPerParts)
    {
    	_abundanceMin = abundanceMin;
    	_abundanceMax = abundanceMax;
    	_bankIndex = bankIndex;
    }

	~SimkaCompressedProcessor(){}
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCompressedProcessor (_bags, _nbKmerPerParts, _nbDistinctKmerPerParts, _chordPerParts, _abundanceMin, _abundanceMax, _bankIndex);  }
    //CountProcessorAbstract<span>* clone ()  {  return new SimkaCompressedProcessor (_bags, _caches, _cacheIndexes, _abundanceMin, _abundanceMax);  }
	void finishClones (vector<ICountProcessor<span>*>& clones){}

	bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum){

		if(count[0] < _abundanceMin || count[0] > _abundanceMax) return false;

		//if(count[0] > 1000)
		//cout << kmer.toString(31) << "  " << count[0] << endl;
		Kmer_BankId_Count item(kmer, _bankIndex, count[0]);
		_bags[partId]->insert(item);
		_nbDistinctKmerPerParts[partId] += 1;
		_nbKmerPerParts[partId] += count[0];
		_chordPerParts[partId] += pow(count[0], 2);

		/*
		size_t index = _cacheIndexes[partId];

		_caches[partId][index] = item;
		index += 1;

		if(index == NB_COUNT_CACHE){
			_bags[partId]->insert(_caches[partId], index);
			_cacheIndexes[partId] = 0;
		}
		else{
			_cacheIndexes[partId] = index;
		}*/

		return true;
	}


	vector<Bag<Kmer_BankId_Count>* >& _bags;
	vector<u_int64_t>& _nbDistinctKmerPerParts;
	vector<u_int64_t>& _nbKmerPerParts;
	vector<u_int64_t>& _chordPerParts;
	CountNumber _abundanceMin;
	CountNumber _abundanceMax;
	size_t _bankIndex;
	//_stats->_chord_N2[i] += pow(abundanceI, 2);
	//vector<vector<Count> >& _caches;
	//vector<size_t>& _cacheIndexes;
};










template<size_t span>
class SimkaAbundanceProcessor : public CountProcessorAbstract<span>{

public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;

    //SimkaCompressedProcessor(vector<BagGzFile<Count>* >& bags, vector<vector<Count> >& caches, vector<size_t>& cacheIndexes, CountNumber abundanceMin, CountNumber abundanceMax) : _bags(bags), _caches(caches), _cacheIndexes(cacheIndexes)
    SimkaAbundanceProcessor(CountNumber abundanceMin, CountNumber abundanceMax)
    {
    	_abundanceMin = abundanceMin;
    	_abundanceMax = abundanceMax;
    	//_bankIndex = bankIndex;
    }

	~SimkaAbundanceProcessor(){}
    CountProcessorAbstract<span>* clone ()  {  return new SimkaAbundanceProcessor<span>(_abundanceMin, _abundanceMax);  }
    //CountProcessorAbstract<span>* clone ()  {  return new SimkaCompressedProcessor (_bags, _caches, _cacheIndexes, _abundanceMin, _abundanceMax);  }
	void finishClones (vector<ICountProcessor<span>*>& clones){}

	bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum){

		if(count[0] < _abundanceMin || count[0] > _abundanceMax) return false;



		/*
		size_t index = _cacheIndexes[partId];
		_caches[partId][index] = item;
		index += 1;
		if(index == NB_COUNT_CACHE){
			_bags[partId]->insert(_caches[partId], index);
			_cacheIndexes[partId] = 0;
		}
		else{
			_cacheIndexes[partId] = index;
		}*/

		return true;
	}

	CountNumber _abundanceMin;
	CountNumber _abundanceMax;
	//_stats->_chord_N2[i] += pow(abundanceI, 2);
	//vector<vector<Count> >& _caches;
	//vector<size_t>& _cacheIndexes;
};


template<size_t span=KMER_DEFAULT_SPAN>
class SimkaPartitionWriter
{
public:


    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    //typedef tuple<Count, StorageItKmerCount<span>*> KmerCount_It;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;


	string _outputDir;
	size_t _nbPartitions;

	vector<u_int64_t> _nbKmerPerParts;
	vector<u_int64_t> _nbDistinctKmerPerParts;
	vector<u_int64_t> _chordNiPerParts;
	vector<Bag<Kmer_BankId_Count>* > _bags;
	vector<Bag<Kmer_BankId_Count>* > _cachedBags;

	SimkaPartitionWriter(const string& oututDir, size_t nbPartitions){
		_outputDir = oututDir;
		_nbPartitions = nbPartitions;

		_nbKmerPerParts = vector<u_int64_t>(_nbPartitions, 0);
		_nbDistinctKmerPerParts =  vector<u_int64_t>(_nbPartitions, 0);
		_chordNiPerParts = vector<u_int64_t>(_nbPartitions, 0);


		//vector<Bag<Kmer_BankId_Count>* > bags;
		//vector<Bag<Kmer_BankId_Count>* > cachedBags;
		for(size_t i=0; i<_nbPartitions; i++){
			//string outputFilename = _outputDir + "/" + _datasetID + "_" + Stringify::format("%i", i) + ".gz";
			string outputFilename = _outputDir + "/" + Stringify::format("%i", i) + ".gz";
			Bag<Kmer_BankId_Count>* bag = new BagGzFile<Kmer_BankId_Count>(outputFilename);
			Bag<Kmer_BankId_Count>* cachedBag = new BagCache<Kmer_BankId_Count>(bag, 10000);
			_cachedBags.push_back(cachedBag);
			//BagCache bagCache(*bag, 10000);
			_bags.push_back(bag);
		}
	}

	void insert(Type& kmer, u_int64_t bankId, u_int64_t abundance){

		size_t part = oahash(kmer) % _nbPartitions;
		_cachedBags[part]->insert(Kmer_BankId_Count(kmer, bankId, abundance));
		_nbDistinctKmerPerParts[part] += 1;
		_nbKmerPerParts[part] += abundance;
		_chordNiPerParts[part] += pow(abundance, 2);

	}

	void end(){
		for(size_t i=0; i<_nbPartitions; i++){
			//bags[i]->flush();
			//cachedBags[i]->flush();
			delete _cachedBags[i];
			//delete bags[i];
		}
	}
};
/*

class SimkaCompressedProcessor_Mini{

public:

    typedef typename Kmer<>::Type  Type;
    typedef typename Kmer<>::Count Count;

    //SimkaCompressedProcessor(vector<BagGzFile<Count>* >& bags, vector<vector<Count> >& caches, vector<size_t>& cacheIndexes, CountNumber abundanceMin, CountNumber abundanceMax) : _bags(bags), _caches(caches), _cacheIndexes(cacheIndexes)
    SimkaCompressedProcessor_Mini(vector<BagGzFile<Count>* >& bags, vector<u_int64_t>& nbKmerPerParts, vector<u_int64_t>& nbDistinctKmerPerParts, vector<u_int64_t>& chordPerParts, CountNumber abundanceMin, CountNumber abundanceMax) :
    	_bags(bags), _nbDistinctKmerPerParts(nbDistinctKmerPerParts), _nbKmerPerParts(nbKmerPerParts), _chordPerParts(chordPerParts)
    {
    	_abundanceMin = abundanceMin;
    	_abundanceMax = abundanceMax;
    }


	bool process (size_t partId, const Type& kmer, CountType count){

		if(count < _abundanceMin || count > _abundanceMax) return false;

		Count item(kmer, count);
		_bags[partId]->insert(item);
		_nbDistinctKmerPerParts[partId] += 1;
		_nbKmerPerParts[partId] += count;
		_chordPerParts[partId] += pow(count, 2);



		return true;
	}


	vector<BagGzFile<Count>* >& _bags;
	vector<u_int64_t>& _nbDistinctKmerPerParts;
	vector<u_int64_t>& _nbKmerPerParts;
	vector<u_int64_t>& _chordPerParts;
	CountNumber _abundanceMin;
	CountNumber _abundanceMax;
	//_stats->_chord_N2[i] += pow(abundanceI, 2);
	//vector<vector<Count> >& _caches;
	//vector<size_t>& _cacheIndexes;
};*/





template<size_t span>
class SimkaMiniKmerCounter : public Algorithm{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    typedef typename Kmer<span>::ModelCanonical                             Model;
    typedef typename Kmer<span>::ModelCanonical::Iterator                             ModelIt;


    //typedef  Kmer<>::ModelCanonical              ModelCanon;
    //typedef  Kmer<>::ModelMinimizer<ModelCanon>  ModelMini;
    typedef typename Kmer<span>::template ModelMinimizer<Model> ModelMinimizer;

    //typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    //typedef typename ModelCanonical::Kmer                                   KmerType;
    /*
    typedef Kmer<span>::Count           Count;
    typedef Kmer<span>::Type            Type;
    typedef Kmer<span>::ModelCanonical              ModelCanon;
    typedef Kmer<span>::ModelMinimizer<ModelCanon>  Model;*/

	IBank* _bank;
	size_t _kmerSize;
	CountVector* _counts;
    u_int64_t _nbReads;
    string _outputDir;
    size_t _nbPartitions;
    u_int64_t _abundanceMin;
    u_int64_t _abundanceMax;
    SimkaPartitionWriter<span>* _partitionWriter;

	SimkaMiniKmerCounter(IProperties* options, size_t kmerSize, IBank* bank, string outputDir, size_t nbPartitions, u_int64_t abundanceMin, u_int64_t abundanceMax, SimkaPartitionWriter<span>* partitionWriter):
		Algorithm("minikc", -1, options)
	{
		_bank = bank;
		_kmerSize = kmerSize;
		_outputDir = outputDir;
		_nbPartitions = nbPartitions;
		_abundanceMin = abundanceMin;
		_abundanceMax = abundanceMax;
		_partitionWriter = partitionWriter;

		u_int64_t nbCounts = pow(4, _kmerSize);
		//cout << "Nb distinct kmers (canonical): " << nbCounts << endl;
		_counts = new CountVector(nbCounts, 0);
	}

	void execute(){

		count();
		dump();

	}

	void count(){

		_nbReads = 0;
		Iterator<Sequence>* itSeq = createIterator(_bank->iterator(), _bank->estimateNbItems(), "Counting");

		//Model definition of a kmer iterator (this one put kmer in cannonical form)
		//ModelCanonical _model(_kmerSize);
		//Model::
		//Model _kmerIt(_model);
		Model model (_kmerSize);

        // We declare an iterator on a given sequence.
		ModelIt _kmerIt (model);

		Sequence* sequence;

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

			_nbReads += 1;

			sequence = &itSeq->item();

			_kmerIt.setData (sequence->getData());

			for (_kmerIt.first(); !_kmerIt.isDone(); _kmerIt.next()){

				//u_int64_t kmer = min(_kmerIt->value(), revcomp(_kmerIt->value(), _kmerSize)).getVal();
				//Kmer<> canonicalkmer = min(_kmerIt.item(), revcomp(_kmerIt->value()));
				//cout << _kmerIt->value().toString(kmerSize) << endl;

				u_int64_t kmer = _kmerIt->value().getVal();
				//cout << _model.toString(kmer) << endl;
				//cout << kmer << endl;
				(*_counts)[kmer] += 1;

				//cout << kmer << endl;
			}
		}

	}

	void dump(){

		//SimkaPartitionWriter<span>* partitionWriter = new SimkaPartitionWriter<span>(_outputDir, _nbPartitions);

		//ModelMinimizer model (_kmerSize, 7);
		Type kmer;

		//Kmer<>::ModelCanonical _model(_kmerSize);
		//CountVector vec(1, 0);

		for(size_t i=0; i<_counts->size(); i++){

			CountNumber count = (*_counts)[i];
			if(count == 0) continue;
			if(count < _abundanceMin || count > _abundanceMax) continue;


			kmer.setVal(i);

			_partitionWriter->insert(kmer, 0, count);
			//cout << i << " " << model.toString(kmer) << endl;
			//Type kmer(i);
            //u_int64_t mini = model.getMinimizerValue(kmer);
			//size_t p = this->_repartition (mini);

			//vec[0] = count;
			//_proc->process(p, kmer, vec, count);
		}

		_partitionWriter->end();
		//delete partitionWriter;
	}


};



#endif /* GATB_SIMKA_SRC_MINIKC_MINIKC_HPP_ */
