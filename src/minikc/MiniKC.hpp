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

    //SimkaCompressedProcessor(vector<BagGzFile<Count>* >& bags, vector<vector<Count> >& caches, vector<size_t>& cacheIndexes, CountNumber abundanceMin, CountNumber abundanceMax) : _bags(bags), _caches(caches), _cacheIndexes(cacheIndexes)
    SimkaCompressedProcessor(vector<u_int64_t>& nbKmerPerParts, vector<u_int64_t>& nbDistinctKmerPerParts, vector<u_int64_t>& chordPerParts, CountNumber abundanceMin, CountNumber abundanceMax) :
    	_nbDistinctKmerPerParts(nbDistinctKmerPerParts), _nbKmerPerParts(nbKmerPerParts), _chordPerParts(chordPerParts)
    {
    	_abundanceMin = abundanceMin;
    	_abundanceMax = abundanceMax;
    }

	~SimkaCompressedProcessor(){}
    CountProcessorAbstract<span>* clone ()  {  return new SimkaCompressedProcessor (_nbKmerPerParts, _nbDistinctKmerPerParts, _chordPerParts, _abundanceMin, _abundanceMax);  }
    //CountProcessorAbstract<span>* clone ()  {  return new SimkaCompressedProcessor (_bags, _caches, _cacheIndexes, _abundanceMin, _abundanceMax);  }
	void finishClones (vector<ICountProcessor<span>*>& clones){}

	bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum){

		if(count[0] < _abundanceMin || count[0] > _abundanceMax) return false;

		//Count item(kmer, count[0]);
		//_bags[partId]->insert(item);
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


	//vector<BagGzFile<Count>* >& _bags;
	vector<u_int64_t>& _nbDistinctKmerPerParts;
	vector<u_int64_t>& _nbKmerPerParts;
	vector<u_int64_t>& _chordPerParts;
	CountNumber _abundanceMin;
	CountNumber _abundanceMax;
	//_stats->_chord_N2[i] += pow(abundanceI, 2);
	//vector<vector<Count> >& _caches;
	//vector<size_t>& _cacheIndexes;
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
class MiniKC : public Algorithm{

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
    Repartitor& _repartition;
    SimkaCompressedProcessor<span>* _proc;
    u_int64_t _nbReads;

	MiniKC(IProperties* options, size_t kmerSize, IBank* bank, Repartitor& repartition, SimkaCompressedProcessor<span>* proc):
		Algorithm("minikc", -1, options), _repartition(repartition)
	{
		_bank = bank;
		_kmerSize = kmerSize;
		_proc = proc;


		u_int64_t nbCounts = pow(4, _kmerSize);
		cout << "Nb distinct kmers (canonical): " << nbCounts << endl;
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

		ModelMinimizer model (_kmerSize, 7);
		Type kmer;

		//Kmer<>::ModelCanonical _model(_kmerSize);
		CountVector vec(1, 0);

		for(size_t i=0; i<_counts->size(); i++){

			CountNumber count = (*_counts)[i];
			if(count == 0) continue;

			kmer.setVal(i);

			//cout << i << " " << model.toString(kmer) << endl;
			//Type kmer(i);
            u_int64_t mini = model.getMinimizerValue(kmer);
			size_t p = this->_repartition (mini);

			vec[0] = count;
			_proc->process(p, kmer, vec, count);
		}

	}


};



#endif /* GATB_SIMKA_SRC_MINIKC_MINIKC_HPP_ */
