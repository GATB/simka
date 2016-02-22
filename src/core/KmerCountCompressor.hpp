
#ifndef _GATB_CORE_KMER_IMPL_KMER_COUNT_COMPRESSOR_HPP_
#define _GATB_CORE_KMER_IMPL_KMER_COUNT_COMPRESSOR_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/compression/RangeCoder.hpp>
#include <gatb/tools/compression/CompressionUtils.hpp>
#include <gatb/system/impl/System.hpp>


#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

const u_int64_t MAX_MEMORY_PER_BLOCK = 100000;
//#define INDEXING
#define MONO_BANK

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/




/*********************************************************************
* ** KmerCountCoder
*********************************************************************/
class KmerCountCoder
{

public:

	KmerCountCoder(int nbBanks, int partitionIndex)
		//_bankCountDeltaModel(3)
	{

		_nbBanks = nbBanks;
		_partitionIndex = partitionIndex;

		_nbKmers = 0;

		/*
    	for(int i=0; i<nbBanks; i++){
    		vector<Order0Model> models;
    		vector<Order0Model> bankModels;
        	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++){
        		models.push_back(Order0Model(256));
        		bankModels.push_back(Order0Model(256));
        	}
    		_abundanceModels.push_back(models);
    		_bankModels.push_back(bankModels);

        	_deltaModels.push_back(Order0Model(3));
        	_bankDeltaModels.push_back(Order0Model(3));
    	}*/

    	//_lastKmerValue = 0;
    	//_lastAbundances.resize(nbBanks, 0);
    	//_lastBanks.resize(nbBanks, 0);


    	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
    		_kmerModel.push_back(Order0Model(256));
    		_bankCountModel.push_back(Order0Model(256));
    	}

		clear();
	}


	void clear(){
    	_lastKmerValue = 0;
    	_lastNbBankCount = 0;

    	//_lastAbundances.clear();
    	//_lastBanks.clear();
		_abundanceModels.clear();
		_bankModels.clear();
		//_deltaModels.clear();
		//_bankDeltaModels.clear();

    	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
    		_kmerModel[i].clear();
    		_bankCountModel[i].clear();
    	}

    	//_bankCountDeltaModel.clear();
	}

	void addField(){
		//_lastAbundances.push_back(0);
		//_lastBanks.push_back(0);

		vector<Order0Model> models;
		vector<Order0Model> bankModels;
    	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++){
    		models.push_back(Order0Model(256));
    		bankModels.push_back(Order0Model(256));
    	}
		_abundanceModels.push_back(models);
		_bankModels.push_back(bankModels);

		//_deltaModels.push_back(Order0Model(3));
		//_bankDeltaModels.push_back(Order0Model(3));


	}

	u_int64_t getNbKmers(){
		return _nbKmers;
	}

    static string toString(u_int64_t value){
    	char buffer[40];
    	snprintf(buffer, 30, "%llu", value);
    	return string(buffer);
    }

protected:

    int _nbBanks;
	int _partitionIndex;
	u_int64_t _nbKmers;

	vector<Order0Model> _kmerModel;
	vector<vector<Order0Model> > _bankModels;
	vector<vector<Order0Model> > _abundanceModels;
	//vector<Order0Model> _bankDeltaModels;
	//vector<Order0Model> _deltaModels;
	//vector<u_int16_t> _lastAbundances;
	//vector<u_int16_t> _lastBanks;

    u_int64_t _lastKmerValue;
    u_int64_t _lastNbBankCount;

    //vector<u_int64_t> _blockSizes;

	vector<Order0Model> _bankCountModel;
	//Order0Model _bankCountDeltaModel;


};



/*********************************************************************
* ** KmerCountCompressorPartition
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class KmerCountCompressorPartition : public KmerCountCoder
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    KmerCountCompressorPartition(const string& outputDir, int partitionIndex, int nbBanks)
    : KmerCountCoder(nbBanks, partitionIndex)
    {
    	string filename = outputDir + "/part_" + KmerCountCoder::toString(_partitionIndex);
    	_outputFile = System::file().newFile(filename.c_str(), "wb");
    }

    ~KmerCountCompressorPartition(){




    	//cout << _rangeEncoder.getBufferSize() << endl;


		/*
		string path = _outputFile->getPath();
		cout << "Partition " << _partitionIndex << endl;
		cout << "\tNb kmers: " << _nbKmers << endl;
		cout << "\tCompressed size: " << System::file().getSize(path) << endl;
		cout << "\tByte per kmer count: " << System::file().getSize(path) / (float)_nbKmers<< endl;*/


		delete _outputFile;
    }

    void flush(){
    	_rangeEncoder.flush();

    	writeBlock();

    	clear();
    	_rangeEncoder.clear();
    	CompressionUtils::encodeNumeric(_rangeEncoder, _kmerModel, _nbKmers);
    	_rangeEncoder.flush();
    	//for(u_int64_t blockSize : _blockSizes){
    	//	CompressionUtils::encodeNumeric(_rangeEncoder, _kmerModel, blockSize);
    	//}
    	_outputFile->fwrite((const char*) _rangeEncoder.getBuffer(true), _rangeEncoder.getBufferSize(), 1);

		_outputFile->flush();
    }

    void insert(const Type& kmer, const CountVector& abundancePerBank){

    	_nbKmers += 1;

    	u_int64_t kmerValue = kmer.getVal();
    	CompressionUtils::encodeNumeric(_rangeEncoder, _kmerModel, kmerValue - _lastKmerValue);
    	_lastKmerValue = kmerValue;


    	if(abundancePerBank.size() == 1){

        	CompressionUtils::encodeNumeric(_rangeEncoder, _bankCountModel, abundancePerBank[0]);
    	}
    	else{





			//u_int64_t deltaValue;
			//u_int8_t deltaType;

			int modelIndex = 0;

			_banks.clear();
			//int nbBankCount;
			for(size_t bankId=0; bankId<abundancePerBank.size(); bankId++){
				if(abundancePerBank[bankId] > 0){
					_banks.push_back(bankId);
					//nbBankCount += 1;
				}
			}

			//deltaType = CompressionUtils::getDeltaValue(nbBankCount, _lastNbBankCount, &deltaValue);
			//_rangeEncoder.encode(_bankCountDeltaModel, deltaType);
			//CompressionUtils::encodeNumeric(_rangeEncoder, _bankCountModel, deltaValue);
			//_lastNbBankCount = nbBankCount;
			CompressionUtils::encodeNumeric(_rangeEncoder, _bankCountModel, _banks.size());

			int lastBankId = 0;

			for(size_t i=0; i<_banks.size(); i++){

				int bankId = _banks[i];

				u_int16_t abundance = abundancePerBank[bankId];

				//if(abundance == 0){

				//}
				//else{


					if(modelIndex >= _bankModels.size()){
						addField();
					}

					//deltaType = CompressionUtils::getDeltaValue(bankId, _lastBanks[modelIndex], &deltaValue);
					//_rangeEncoder.encode(_bankDeltaModels[modelIndex], deltaType);
					//CompressionUtils::encodeNumeric(_rangeEncoder, _bankModels[modelIndex], deltaValue);
					//_lastBanks[modelIndex] = bankId;

					CompressionUtils::encodeNumeric(_rangeEncoder, _bankModels[modelIndex], bankId - lastBankId);
					lastBankId = bankId;

					//deltaType = CompressionUtils::getDeltaValue(abundance, _lastAbundances[modelIndex], &deltaValue);
					//_rangeEncoder.encode(_deltaModels[modelIndex], deltaType);
					//CompressionUtils::encodeNumeric(_rangeEncoder, _abundanceModels[modelIndex], deltaValue);
					//_lastAbundances[modelIndex] = abundance;
					CompressionUtils::encodeNumeric(_rangeEncoder, _abundanceModels[modelIndex], abundance);

					modelIndex += 1;

					//}

			}
    	}


    	if(_rangeEncoder.getBufferSize() >= MAX_MEMORY_PER_BLOCK){
    		writeBlock();
    	}

    }


	void writeBlock(){
    	if(_rangeEncoder.getBufferSize() > 0){
    		//_rangeEncoder.flush();
			//_blockSizes.push_back(_rangeEncoder.getBufferSize());
			_outputFile->fwrite((const char*) _rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), 1);
    	}
    	_rangeEncoder.clearBuffer();
		//_rangeEncoder.clear();
		//clear();
	}

	u_int64_t getSizeByte(){
		//cout << _outputFile->getPath() << endl;
		//cout << System::file().getSize(_outputFile->getPath()) << endl;
		return System::file().getSize(_outputFile->getPath());
	}



private:

    RangeEncoder _rangeEncoder;
    IFile* _outputFile;
    vector<int> _banks;

};





/*********************************************************************
* ** KmerCountCompressor
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class KmerCountCompressor
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    /** */
    KmerCountCompressor(const string& outputDir, int nbPartitions, int nbBanks){

    	_rootDir = outputDir; // + "/dsk_output/";
    	_nbPartitions = nbPartitions;
    	_nbBanks = nbBanks;

		System::file().rmdir(_rootDir);
		System::file().mkdir(_rootDir, -1);

		cout << _rootDir << endl;
		//cout << nbPartitions << endl;
		for(int i=0; i<nbPartitions; i++){
			KmerCountCompressorPartition<span>* comp = new KmerCountCompressorPartition<span>(_rootDir, i, nbBanks);
			_partitionCompressors.push_back(comp);
		}
    }

    ~KmerCountCompressor(){

    	u_int64_t nbKmers = 0;
    	u_int64_t size = 0;

    	for(size_t i=0; i<_partitionCompressors.size(); i++){
    		KmerCountCompressorPartition<span>* comp = _partitionCompressors[i];
    		comp->flush();
    		nbKmers += comp->getNbKmers();
    		size += comp->getSizeByte();
    		delete comp;
    	}


		cout << "Compression statistics " << endl;
		cout << "\tNb kmers: " << nbKmers << endl;
		cout << "\tCompressed size: " << size << "B  -  " << size/MBYTE << " MB" << endl;
		cout << "\tByte per kmer count: " << size / (float) nbKmers<< endl;

        IFile* outputFile = System::file().newFile(_rootDir + "/dsk_count_data", "wb");
        outputFile->print("%i %i", _nbPartitions, _nbBanks);
        outputFile->flush();
        delete outputFile;
    }


    void insert(int partitionIndex, const Type& kmer, const CountVector& abundancePerBank){
    	_partitionCompressors[partitionIndex]->insert(kmer, abundancePerBank);
    }


private:

    string _rootDir;
    int _nbPartitions;
    int _nbBanks;

	vector<KmerCountCompressorPartition<span>* > _partitionCompressors;

};






























/*********************************************************************
* ** KmerCountDecompressorPartition
*********************************************************************/
template<typename Functor, size_t span=KMER_DEFAULT_SPAN>
class KmerCountDecompressorPartition : KmerCountCoder
{
public:

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    KmerCountDecompressorPartition(const string& inputDir, int partitionIndex, int nbBanks, Functor* functor, gatb::core::tools::dp::IteratorListener* progress)
    : KmerCountCoder(nbBanks, partitionIndex)
    {

    	_progress = progress;
    	_functor = functor;

    	_nbDecodedKmers = 0;
    	_nbDecodedKmersProgress = 0;

    	string filename = inputDir + "/part_" + KmerCountCoder::toString(_partitionIndex);
    	_inputFile = new ifstream(filename.c_str(), ios::in|ios::binary);

    	_inputFile->seekg(0, _inputFile->end);
    	_rangeDecoder.setInputFile(_inputFile, true);
    	_nbKmers = CompressionUtils::decodeNumeric(_rangeDecoder, _kmerModel);
    	//cout << _nbKmers << endl;
    	//_rangeEncoder.flush();
    	//for(u_int64_t blockSize : _blockSizes){
    	//	CompressionUtils::encodeNumeric(_rangeEncoder, _kmerModel, blockSize);
    	//}
    	//_outputFile->fwrite((const char*) _rangeEncoder.getBuffer(true), _rangeEncoder.getBufferSize(), 1);


    	clear();
    	_rangeDecoder.clear();
    	_inputFile->seekg(0, _inputFile->beg);
    	_rangeDecoder.setInputFile(_inputFile);


    	/*
    	for(int i=0; i<nbBanks; i++){
    		vector<Order0Model> models;
    		vector<Order0Model> bankModels;
        	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++){
        		models.push_back(Order0Model(256));
        		bankModels.push_back(Order0Model(256));
        	}
    		_abundanceModels.push_back(models);
    		_bankModels.push_back(bankModels);

        	_deltaModels.push_back(Order0Model(3));
        	_bankDeltaModels.push_back(Order0Model(3));
    	}

    	_lastKmerValue = 0;
    	_lastAbundances.resize(nbBanks, 0);
    	_lastBanks.resize(nbBanks, 0);

    	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
    		_kmerModel.push_back(Order0Model(256));
    	}*/
    }

    ~KmerCountDecompressorPartition(){

    	//delete _functor;
    	delete _inputFile;

    }

    void execute(){

    	while(_nbDecodedKmers < _nbKmers){

    		//cout << _nbDecodedKmers << endl;

        	u_int64_t kmerDeltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _kmerModel);
        	u_int64_t kmerValue = _lastKmerValue + kmerDeltaValue;

        	_lastKmerValue = kmerValue;

        	CountVector abundancePerBanks;

        	u_int64_t nbCounts = CompressionUtils::decodeNumeric(_rangeDecoder, _bankCountModel);
        	while(_bankModels.size() < nbCounts){
        		addField();
        	}

        	u_int64_t lastBankId = 0;
			u_int64_t bankId;

			for(int modelIndex=0; modelIndex<nbCounts; modelIndex++){

				u_int64_t bankIdDelta = CompressionUtils::decodeNumeric(_rangeDecoder, _bankModels[modelIndex]);
				bankId = lastBankId + bankIdDelta;
				lastBankId = bankId;

				while(abundancePerBanks.size() < bankId){
					abundancePerBanks.push_back(0);
				}

				//deltaType = CompressionUtils::getDeltaValue(abundance, _lastAbundances[modelIndex], &deltaValue);
				//_rangeEncoder.encode(_deltaModels[modelIndex], deltaType);
				//CompressionUtils::encodeNumeric(_rangeEncoder, _abundanceModels[modelIndex], deltaValue);
				//_lastAbundances[modelIndex] = abundance;
				u_int64_t abundance = CompressionUtils::decodeNumeric(_rangeDecoder, _abundanceModels[modelIndex]);
				abundancePerBanks.push_back(abundance);

			}

			while(abundancePerBanks.size() < _nbBanks){
				abundancePerBanks.push_back(0);
			}

			Type kmer(kmerValue);
	    	_functor->execute(kmer, abundancePerBanks);

	    	_nbDecodedKmers += 1;
	    	_nbDecodedKmersProgress += 1;

	        if (_nbDecodedKmersProgress > 500000)   {  _progress->inc (_nbDecodedKmersProgress); _nbDecodedKmersProgress = 0;  }
    	}

    	_progress->inc (_nbDecodedKmersProgress);
    }


private:

    Functor* _functor;
    gatb::core::tools::dp::IteratorListener* _progress;

    RangeDecoder _rangeDecoder;
    ifstream* _inputFile;

    u_int64_t _nbDecodedKmers;
    u_int64_t _nbDecodedKmersProgress;
};


/*********************************************************************
* ** KmerCountDecompressor
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class KmerCountDecompressor : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    /** */
    KmerCountDecompressor(const string& inputDir, int nbCores)
    : Algorithm("kcc", 0, 0)
    {

    	//getInput()->setStr(STR_VERBOSE, "1");

    	_inputDir = inputDir;
    	_nbCores = nbCores;
    	IFile* dskCountDataFile = System::file().newFile(inputDir + "/dsk_count_data", "rb");

    	vector<string> numbers;
		string n = "";
    	while(true){

    		u_int8_t c = dskCountDataFile->get();
    		if(c == ' ' || dskCountDataFile->isEOF()){
    			numbers.push_back(n);
    			n.clear();
    		}
    		else{
    			n += c;
    		}

    		if(dskCountDataFile->isEOF()) break;
    	}

    	//for(string& number: numbers){
    	//	cout << number << endl;
    	//}

    	_nbPartitions = atoi(numbers[0].c_str());
    	_nbBanks = atoi(numbers[1].c_str());

    	//for(int i=0; i<_nbCores; i++)

		//for(int i=0; i<_nbPartitions; i++){
		//	KmerCountCompressorPartition<span> decomp(inputDir, i, _nbBanks);
		//	decomp.iterKmers(_nbCores);
		//}

    }


    ~KmerCountDecompressor(){
        delete _progress;
    }

    void setupProgress(){

        RangeDecoder rangeDecoder;
    	vector<Order0Model> kmerModel;
    	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
    		kmerModel.push_back(Order0Model(256));
    	}

    	u_int64_t totalKmers = 0;

    	for(int i=0; i<_nbPartitions; i++){

    		string filename = _inputDir + "/part_" + KmerCountCoder::toString(i);
    		ifstream* file = new ifstream(filename.c_str(), ios::in|ios::binary);

    		file->seekg(0, file->end);
    		rangeDecoder.setInputFile(file, true);
        	u_int64_t nbKmers = CompressionUtils::decodeNumeric(rangeDecoder, kmerModel);
        	totalKmers += nbKmers;

        	rangeDecoder.clear();
        	for(int j=0; j<CompressionUtils::NB_MODELS_PER_NUMERIC; j++){
        		kmerModel[j].clear();
        	}
    		delete file;
    	}

    	_progress = new gatb::core::tools::misc::impl::ProgressSynchro (
            createIteratorListener (totalKmers, "Decompressing kmer counts"),
            System::thread().newSynchronizer());
        _progress->init ();
    }

    template <typename Functor>
    static void *callMyFunction(void *object){
		((KmerCountDecompressorPartition<Functor, span>*)object)->execute();
		//object->execute();
		return NULL;
    }

    void execute(){

    }

    template <typename Functor>
    void iterate (const Functor& functor, size_t groupSize=1000){

    	setupProgress();

    	pthread_t* tab_threads = new pthread_t[_nbCores];
    	//thread_arg_decoder *  targ = new thread_arg_decoder [_nbCores];

		vector<KmerCountDecompressorPartition<Functor, span>* > _partitionDecompressors;

    	for(int i=0; i<_nbPartitions;){

    		for(int j=0; j<_nbCores && i<_nbPartitions; j++){
    			Functor* func = new Functor(functor);

    			KmerCountDecompressorPartition<Functor, span>* decomp = new KmerCountDecompressorPartition<Functor, span>(_inputDir, i, _nbBanks, func, _progress);

    			_partitionDecompressors.push_back(decomp);

    			i += 1;
    		}

    		//Lala * lala = new Lala();
    		for(int j=0; j<_partitionDecompressors.size(); j++){
        		//cout << "start" << endl;
                pthread_create(&tab_threads[j], NULL, &KmerCountDecompressor::callMyFunction<Functor>, _partitionDecompressors[j]);
    			//_partitionDecompressors[j]->execute();
        		//cout << "loulou" << endl;
			}

    		for(int j=0; j<_partitionDecompressors.size(); j++){
    			pthread_join(tab_threads[j], NULL);
			}

    		for(int j=0; j<_partitionDecompressors.size(); j++){
    			delete _partitionDecompressors[j];
    		}

    		_partitionDecompressors.clear();

    	}

    	delete tab_threads;

        _progress->finish ();
    }



private:

    string _inputDir;
    int _nbCores;
    int _nbBanks;
    int _nbPartitions;

    gatb::core::tools::dp::IteratorListener* _progress;
	//vector<KmerCountDecompressorPartition<span>* > _partitionDecompressors;

};





/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_KMER_COUNT_COMPRESSOR__HPP_ */
