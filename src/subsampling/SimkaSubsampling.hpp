/*
 * Subsampling.hpp
 *
 *  Created on: 26 oct. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_SUBSAMPLING_SUBSAMPLING_HPP_
#define GATB_SIMKA_SRC_SUBSAMPLING_SUBSAMPLING_HPP_

#include <gatb/gatb_core.hpp>




class SimkaSubsampling {

public:

	size_t _nbBanks;
	string _outputDir;
	vector<string> _bankNames;
	size_t _kmerSize;




	SimkaSubsampling(const string& outputDir, size_t nbBanks, vector<string>& bankNames, size_t kmerSize):
		_bankNames(bankNames)
	{
		_outputDir = outputDir;
		_nbBanks = nbBanks;
		_kmerSize = kmerSize;
	}

	SimkaSubsampling(const string& outputDir, size_t kmerSize)
	{
		_outputDir = outputDir;
		_kmerSize = kmerSize;
	}

	~SimkaSubsampling(){

	}

	typedef tuple<u_int64_t, size_t> BankSize_BankId;
	static bool sortFileBySize (BankSize_BankId& i, BankSize_BankId& j){
		return ( get<0>(i) < get<0>(j) );
	}

	void setup(){

		string inputDir = _outputDir + "/input/";

		vector<BankSize_BankId> filenameSizes;

		for(size_t i=0; i<_nbBanks; i++){

			IBank* bank = Bank::open(inputDir + _bankNames[i]);
			LOCAL(bank);
			u_int64_t nbReads;
			u_int64_t readSize;
			u_int64_t totalSize;
			bank->estimate (nbReads, totalSize, readSize);

			//cout << i << "    " << nbReads << " " << totalSize << " " << readSize << endl;
			filenameSizes.push_back(BankSize_BankId(totalSize, i));
		}

		sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);

		size_t smallerBankId = get<1>(filenameSizes[0]);

		IBank* bank = Bank::open(inputDir + _bankNames[smallerBankId]);
		LOCAL(bank);

		u_int64_t nbKmers = 0;
		u_int64_t realSize = 0;
		Iterator<Sequence>* it = bank->iterator();
		for(it->first(); !it->isDone(); it->next()){
			if(it->item().getDataSize() < _kmerSize) continue;

			nbKmers += it->item().getDataSize() - _kmerSize + 1;
			//cout << it->item().getDataSize() << endl;
		}

		cout << "Subsampling max: " << nbKmers << endl;
		//cout << "Subsampling Space Size: " << nbKmers << endl;
		//cout << smallerBankId << " " << realSize << endl;
	}

	void setupRandomNumberGenerator(u_int64_t min, u_int64_t max){
		srand( time( NULL ) );
		_rngMin = min;
		_rngN = (int)(max - min + 1);
		_rngRemainder = RAND_MAX % _rngN;
		//int x, output;
	}

	int generateUniformInt(){
		int x;
		int output;
		do {
		  x = rand();
		  output = x % _rngN;
		} while (x >= RAND_MAX - _rngRemainder);
		return _rngMin + output;
	}


	u_int64_t _rngRemainder;
	u_int64_t _rngMin;
	int _rngN;

	void sampleWithReplacement(const string& bankId, u_int64_t maxKmersToPick, u_int64_t maxKmersSpace, vector<u_int32_t>& pickedReads){

		vector<u_int32_t> nbKmersPerReads;

		string inputDir = _outputDir + "/input/";
		IBank* bank = Bank::open(inputDir + bankId);
		LOCAL(bank);

		u_int64_t nbKmers = 0;
		u_int64_t realSize = 0;
		Iterator<Sequence>* it = bank->iterator();
		for(it->first(); !it->isDone(); it->next()){
			if(it->item().getDataSize() < _kmerSize){
				nbKmersPerReads.push_back(0);
			}
			else{
				u_int64_t nbKmersInRead = it->item().getDataSize() - _kmerSize + 1;

				nbKmersPerReads.push_back(nbKmersInRead);
				nbKmers += nbKmersInRead;
			}

			if(nbKmers > maxKmersSpace){
				break;
			}
		}

		//cout << "subsampling" << endl;
		//cout << nbKmersPerReads.size() << endl;
		setupRandomNumberGenerator(0, nbKmersPerReads.size());
		//cout << _nbKmersPerReads.size() << endl;

		nbKmers = 0;
		//vector<u_int32_t> pickedReads(nbKmersPerReads.size(), 0);
		pickedReads.clear();
		pickedReads.resize(nbKmersPerReads.size(), 0);

		u_int64_t nbPickedReads = 0;

		//cout << "max to pick:   " << maxKmersToPick << endl;

		while(true){

			//cout << nbKmers << "   " << maxKmersToPick << endl;
			//if(nbPickedReads < maxReads){
			//	break;
			//}
			if(nbKmers > maxKmersToPick){
				break;
			}

			int r = generateUniformInt();
			pickedReads[r] += 1;
			nbKmers += nbKmersPerReads[r];

		}

		//return pickedReads;
		//for(size_t i=0; i<pickedReads.size(); i++){
		//	if(pickedReads[i] > 0){
		//		cout << i << "   " << pickedReads[i] << endl;
		//	}
		//}
	}



};

#endif /* GATB_SIMKA_SRC_SUBSAMPLING_SUBSAMPLING_HPP_ */
