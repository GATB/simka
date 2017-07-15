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
	u_int64_t _nbKmersToPick;




	SimkaSubsampling(const string& outputDir, size_t nbBanks, vector<string>& bankNames, size_t kmerSize):
		_bankNames(bankNames)
	{
		_outputDir = outputDir;
		_nbBanks = nbBanks;
		_kmerSize = kmerSize;
		_nbKmersToPick = 0;
	}

	SimkaSubsampling(const string& outputDir, size_t kmerSize)
	{
		_outputDir = outputDir;
		_kmerSize = kmerSize;
		_nbKmersToPick = 0;
	}

	~SimkaSubsampling(){

	}

	typedef tuple<u_int64_t, size_t> BankSize_BankId;
	//static bool sortFileBySize (BankSize_BankId& i, BankSize_BankId& j){
	//	return ( get<0>(i) < get<0>(j) );
	//}
	struct SortFileBySize
	{
	    bool operator()(const BankSize_BankId& i, const BankSize_BankId& j) const
	    {
	    	return ( get<0>(i) < get<0>(j) );
	    }
	};

	void setup(){

		string dataset_info_filename = _outputDir + "/" + "subsampling_info.csv";
		ofstream dataset_info_file(dataset_info_filename.c_str());
		//cout << dataset_info_filename << endl;

		string line = "DATASET_ID;NB_READS;MEAN_READ_LENGTH\n";
		dataset_info_file.write(line.c_str(), line.size());

		string inputDir = _outputDir + "/input/";

		u_int64_t totalNbReads = 0;
		u_int64_t totalLengthOfReads = 0;
		u_int64_t minReads = -1;
		u_int64_t maxReads = 0;
		string smallerDatasetId = "";
		string largerDatasetId = "";

		vector<BankSize_BankId> filenameSizes;
		vector<u_int64_t> readLengths;

		for(size_t i=0; i<_nbBanks; i++){

			IBank* bank = Bank::open(inputDir + _bankNames[i]);
			LOCAL(bank);
			u_int64_t nbReads;
			u_int64_t readSize;
			u_int64_t totalSize;
			bank->estimate (nbReads, totalSize, readSize);

			cout << i << ": " << readSize << endl;
			totalNbReads += nbReads;
			totalLengthOfReads += readSize;

			if(nbReads < minReads){
				minReads = nbReads;
				smallerDatasetId = _bankNames[i];
			}
			if(nbReads > maxReads){
				maxReads = nbReads;
				largerDatasetId = _bankNames[i];
			}

			string line = _bankNames[i] + ";" + Stringify::format("%i", nbReads) + ";" + Stringify::format("%i", readSize)  + "\n";
			dataset_info_file.write(line.c_str(), line.size());
			//cout << i << "    " << nbReads << " " << totalSize << " " << readSize << endl;
			filenameSizes.push_back(BankSize_BankId(totalSize, i));
			readLengths.push_back(readSize);
		}
		dataset_info_file.close();

		sort(filenameSizes.begin(),filenameSizes.end(),SortFileBySize());

		//We llok for the dataset which have the closer mean read length than the average read length of all datasets
		u_int64_t minDatasetID = -1;
		int64_t minDistance = -1;
		int64_t meanReadLength = totalLengthOfReads / _nbBanks;
		for(size_t i=0; i<_nbBanks; i++){
			int64_t readLength = readLengths[i];
			u_int64_t distance = abs(readLength - meanReadLength);
			//cout << distance << endl;
			if(distance < minDistance){
				minDistance = distance;
				minDatasetID = i;
			}
		}
		//cout << meanReadLength << endl;


		size_t smallerBankId = get<1>(filenameSizes[0]);

		IBank* bank = Bank::open(inputDir + _bankNames[smallerBankId]);
		LOCAL(bank);

		//u_int64_t nbKmers = 0;
		u_int64_t nbReads = 0;
		Iterator<Sequence>* it = bank->iterator();
		for(it->first(); !it->isDone(); it->next()){
			//if(it->item().getDataSize() < _kmerSize) continue;

			//nbKmers += it->item().getDataSize() - _kmerSize + 1;
			nbReads += 1;
			//cout << it->item().getDataSize() << endl;
		}

		cout << endl << "--------" << endl;
		cout << "Statistics on number of reads in project (estimation): " << endl;
		cout << "\tMean (read per dataset): " << (totalNbReads/_nbBanks) << endl;
		cout << "\tMin (smaller dataset)  : " << minReads  << " (" << smallerDatasetId << ")" << endl;
		cout << "\tMax (larger dataset)   : " << maxReads  << " (" << largerDatasetId << ")" << endl;
		cout << "Complete info filename: " << dataset_info_filename << endl;
		cout << "--------" << endl;

		cout << endl;
		cout << "Reference dataset ID     : " << minDatasetID << endl;
		cout << "Subsampling space (reads): " << nbReads << endl;
		//cout << "Subsampling max: " << nbKmers << endl;
		cout << endl;
		//cout << "Subsampling Space Size: " << nbKmers << endl;
		//cout << smallerBankId << " " << realSize << endl;


		//string setup_filename = _outputDir + "/" + "subsampling_setup.txt";
		//ofstream setup_file(setup_filename.c_str(), ios::binary);
		//dataset_info_file.write((const char*)&nbReads, sizeof(nbReads));
		//dataset_info_file.write((const char*)&minDatasetID, sizeof(minDatasetID));
		//setup_file.close();
	}

	//Compute the maximum number of kmers to pick for the subsample passes
	void start(bool init, u_int64_t referenceDatasetID, u_int64_t nbReadsToPick, u_int64_t subsamplingSpace_maxReads, size_t subsamplingKind){

		//_maxPickableKmers = 0;

		string start_filename = _outputDir + "/" + "subsampling_start.txt";

		if(init){
			_nbKmersToPick = 0;
			//vector<u_int32_t> pickedReads;
			u_int64_t nbPickedKmers = determineSubsamplingSpace(_bankNames[referenceDatasetID]);

			ofstream dataset_info_file(start_filename.c_str(), ios::binary);

			dataset_info_file.write((const char*)&nbPickedKmers, sizeof(nbPickedKmers));
			//string line = "max_kmers:" + Stringify::format("%i", nbPickedKmers);

			dataset_info_file.close();
		}
		else{
			ifstream dataset_info_file(start_filename.c_str(), ios::binary);
			dataset_info_file.read((char*)&_nbKmersToPick, sizeof(_nbKmersToPick));
			dataset_info_file.close();
		}

	}

	u_int64_t determineSubsamplingSpace(const string& bankId){

		//vector<u_int32_t> nbKmersPerReads;

		u_int64_t nbKmers = 0;

		string inputDir = _outputDir + "/input/";
		IBank* bank = Bank::open(inputDir + bankId);
		LOCAL(bank);

		//u_int64_t nbReads = 0;
		//u_int64_t realSize = 0;
		Iterator<Sequence>* it = bank->iterator();
		for(it->first(); !it->isDone(); it->next()){
			if(it->item().getDataSize() < _kmerSize){
				//nbKmersPerReads.push_back(0);
			}
			else{
				u_int64_t nbKmersInRead = it->item().getDataSize() - _kmerSize + 1;
				nbKmers += nbKmersInRead;
				//cout << nbKmersInRead << endl;
				//nbKmersPerReads.push_back(nbKmersInRead);
			}


			//if(nbReads >= subsamplingSpace_maxReads){
			//	break;
			//}

			//nbReads += 1;
		}

		/*
		cout << nbKmersPerReads.size() << endl;

		vector<u_int32_t> readIds(nbKmersPerReads.size(), 0);
		for(size_t i=0; i<readIds.size(); i++){
			readIds[i] = i;
		}
		//cout << readIds.size() << endl;

		//cout << readIds.size() << endl;
		pickedReads.resize(nbKmersPerReads.size(), 0);

		//cout << "subsampling" << endl;
		//cout << nbKmersPerReads.size() << endl;
		setupRandomNumberGenerator(0, nbKmersPerReads.size()-1);
		std::random_shuffle (readIds.begin(), readIds.end(), myrandom);
		//cout << _nbKmersPerReads.size() << endl;

		u_int64_t nbPickedKmers = 0;
		//vector<u_int32_t> pickedReads(nbKmersPerReads.size(), 0);
		pickedReads.clear();
		pickedReads.resize(nbKmersPerReads.size(), 0);

		u_int64_t nbPickedReads = 0;

		//cout << "max to pick:   " << maxKmersToPick << endl;

		for(size_t i=0; i<readIds.size(); i++){

			//cout << nbKmers << "   " << maxKmersToPick << endl;
			//if(nbPickedReads < maxReads){
			//	break;
			//}

			int r = readIds[i];
			//int r = generateUniformInt();
			if(nbKmersPerReads[r] == 0){
				continue;
			}

			pickedReads[r] = 1;
			nbPickedKmers += nbKmersPerReads[r];
			nbPickedReads += 1;

			if(_nbKmersToPick != 0){
				if(nbPickedKmers >= _nbKmersToPick){
					break;
				}
			}
			else if(nbPickedReads >= nbReadsToPick){
				break;
			}

			//cout << r << "   " << nbKmers  << "    " << nbKmersPerReads[r] << endl;
		}

		cout << "nb picked reads: " << nbPickedReads << endl;
		cout << "nb picked k-mers: " << nbPickedKmers << endl;

		//for(size_t i=0; i<pickedReads.size(); i++){
		//	if(pickedReads[i] > 0)
		//	cout << pickedReads[i];
		//}
		//cout << endl;
		 */
		return nbKmers;


	}

	void setupRandomNumberGenerator(u_int64_t min, u_int64_t max){
		srand( time( NULL ) );
		//std::srand ( unsigned ( std::time(0) ) );
		//_rngMin = min;
		_rngN = (int)(max - min + 1);
		//_rngRemainder = RAND_MAX % _rngN;
		//int x, output;
	}

	/*
	int generateUniformInt(){
		int x;
		int output;
		do {
		  x = rand();
		  output = x % _rngN;
		} while (x >= RAND_MAX - _rngRemainder);
		return _rngMin + output;
	}


	int _rngRemainder;
	u_int64_t _rngMin;*/
	int _rngN;

	u_int64_t sample(const string& bankId, u_int64_t nbReadsToPick, u_int64_t subsamplingSpace_maxReads, vector<u_int32_t>& pickedReads, size_t subsamplingKind){
		bool isRandom = (nbReadsToPick != 0);
		cout << "NB read to pick: " << nbReadsToPick << endl;
		cout << "MAX pickable reads: " << subsamplingSpace_maxReads << endl;
		if(subsamplingKind == 0){
			return sampleWithReplacement(bankId, nbReadsToPick, subsamplingSpace_maxReads, pickedReads, isRandom);
		}
		else{
			return sampleWithoutReplacement(bankId, nbReadsToPick, subsamplingSpace_maxReads, pickedReads, isRandom);
		}
	}

	u_int64_t sampleWithReplacement(const string& bankId, u_int64_t nbReadsToPick, u_int64_t subsamplingSpace_maxReads, vector<u_int32_t>& pickedReads, bool isRandom){

		vector<u_int32_t> nbKmersPerReads;

		string inputDir = _outputDir + "/input/";
		IBank* bank = Bank::open(inputDir + bankId);
		LOCAL(bank);

		u_int64_t nbReads = 0;
		u_int64_t realSize = 0;
		Iterator<Sequence>* it = bank->iterator();
		for(it->first(); !it->isDone(); it->next()){
			if(it->item().getDataSize() < _kmerSize){
				nbKmersPerReads.push_back(0);
			}
			else{
				u_int64_t nbKmersInRead = it->item().getDataSize() - _kmerSize + 1;

				//cout << nbKmersInRead << endl;
				nbKmersPerReads.push_back(nbKmersInRead);
			}


			if(nbReads >= subsamplingSpace_maxReads){
				break;
			}

			nbReads += 1;
		}

		cout << nbKmersPerReads.size() << endl;

		//cout << "subsampling" << endl;
		//cout << nbKmersPerReads.size() << endl;
		setupRandomNumberGenerator(0, nbKmersPerReads.size()-1);
		//cout << _nbKmersPerReads.size() << endl;

		u_int64_t nbPickedKmers = 0;
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


			int r = myrandom(_rngN);
			if(nbKmersPerReads[r] == 0){
				continue;
			}

			pickedReads[r] += 1;
			nbPickedKmers += nbKmersPerReads[r];
			nbPickedReads += 1;

			if(_nbKmersToPick != 0){
				if(nbPickedKmers >= _nbKmersToPick){
					break;
				}
			}
			else if(nbPickedReads >= nbReadsToPick){
				break;
			}

			//cout << r << "   " << nbKmers  << "    " << nbKmersPerReads[r] << endl;
		}

		cout << "nb picked reads: " << nbPickedReads << endl;
		cout << "nb picked k-mers: " << nbPickedKmers << endl;


		//for(size_t i=0; i<pickedReads.size(); i++){
		//	if(pickedReads[i] > 0)
		//	cout << pickedReads[i];
		//}
		//cout << endl;

		return nbPickedKmers;
		//if(_nbKmersToPick == 0){
		//	_nbKmersToPick = nbPickedKmers;
		//}
		//cout << endl;

		//return pickedReads;
		//for(size_t i=0; i<pickedReads.size(); i++){
		//	if(pickedReads[i] > 0){
		//		cout << i << "   " << pickedReads[i] << endl;
		//	}
		//}


	}

	static int myrandom (int i) { return std::rand()%i;}

	u_int64_t sampleWithoutReplacement(const string& bankId, u_int64_t nbReadsToPick, u_int64_t subsamplingSpace_maxReads, vector<u_int32_t>& pickedReads, bool isRandom){

		if(isRandom){
			vector<u_int32_t> nbKmersPerReads;

			string inputDir = _outputDir + "/input/";
			IBank* bank = Bank::open(inputDir + bankId);
			LOCAL(bank);

			u_int64_t nbPickedKmers = 0;
			//u_int64_t realSize = 0;
			Iterator<Sequence>* it = bank->iterator();
			for(it->first(); !it->isDone(); it->next()){
				if(it->item().getDataSize() < _kmerSize){
					nbKmersPerReads.push_back(0);
				}
				else{
					u_int64_t nbKmersInRead = it->item().getDataSize() - _kmerSize + 1;
					nbPickedKmers += nbKmersInRead;
					//cout << nbKmersInRead << endl;
					nbKmersPerReads.push_back(nbKmersInRead);
				}


				if(nbPickedKmers >= subsamplingSpace_maxReads){
					break;
				}

				//nbReads += 1;
			}

			//cout << nbKmersPerReads.size() << endl;

			vector<u_int32_t> readIds(nbKmersPerReads.size(), 0);
			for(size_t i=0; i<readIds.size(); i++){
				readIds[i] = i;
			}
			//cout << readIds.size() << endl;

			//cout << readIds.size() << endl;
			//pickedReads.resize(nbKmersPerReads.size(), 0);

			//cout << "subsampling" << endl;
			//cout << nbKmersPerReads.size() << endl;
			setupRandomNumberGenerator(0, nbKmersPerReads.size()-1);
			std::random_shuffle (readIds.begin(), readIds.end(), myrandom);
			//cout << _nbKmersPerReads.size() << endl;

			nbPickedKmers = 0;
			//vector<u_int32_t> pickedReads(nbKmersPerReads.size(), 0);
			pickedReads.clear();
			pickedReads.resize(nbKmersPerReads.size(), 0);

			//u_int64_t nbPickedReads = 0;

			//cout << "max to pick:   " << maxKmersToPick << endl;

			for(size_t i=0; i<readIds.size(); i++){

				//cout << nbKmers << "   " << maxKmersToPick << endl;
				//if(nbPickedReads < maxReads){
				//	break;
				//}

				int r = readIds[i];
				//int r = generateUniformInt();
				if(nbKmersPerReads[r] == 0){
					continue;
				}

				pickedReads[r] = 1;
				nbPickedKmers += nbKmersPerReads[r];
				//nbPickedReads += 1;


				if(nbPickedKmers >= nbReadsToPick){
					break;
				}

				//if(_nbKmersToPick != 0){
				//	if(nbPickedKmers >= _nbKmersToPick){
				//		break;
				//	}
				//}
				//else if(nbPickedReads >= nbReadsToPick){
				//	break;
				//}

				//cout << r << "   " << nbKmers  << "    " << nbKmersPerReads[r] << endl;
			}

			//cout << "nb picked reads: " << nbPickedReads << endl;
			cout << "nb picked k-mers: " << nbPickedKmers << endl;

			//for(size_t i=0; i<pickedReads.size(); i++){
			//	if(pickedReads[i] > 0)
			//	cout << pickedReads[i];
			//}
			//cout << endl;

			return nbPickedKmers;
		}
		else{

			u_int64_t nbPickedKmers = 0;
			//vector<u_int32_t> nbKmersPerReads;

			string inputDir = _outputDir + "/input/";
			IBank* bank = Bank::open(inputDir + bankId);
			LOCAL(bank);

			//u_int64_t readIndex = 0;
			//u_int64_t realSize = 0;
			Iterator<Sequence>* it = bank->iterator();
			for(it->first(); !it->isDone(); it->next()){
				if(it->item().getDataSize() < _kmerSize){
					pickedReads.push_back(0);
				}
				else{
					u_int64_t nbKmersInRead = it->item().getDataSize() - _kmerSize + 1;
					nbPickedKmers += nbKmersInRead;
					pickedReads.push_back(1);

					//cout << nbKmersInRead << endl;
					//nbKmersPerReads.push_back(nbKmersInRead);
				}


				if(nbPickedKmers >= subsamplingSpace_maxReads){
					break;
				}

				//readIndex += 1;
			}


			return nbPickedKmers;
		}

		//if(_nbKmersToPick == 0){
		//	_nbKmersToPick = nbPickedKmers;
		//}
		//cout << endl;
		//for(size_t i=0; i<pickedReads.size(); i++){
		//	cout << pickedReads[i];
		//}
		//cout << endl;

		//return pickedReads;
		//for(size_t i=0; i<pickedReads.size(); i++){
		//	if(pickedReads[i] > 0){
		//		cout << i << "   " << pickedReads[i] << endl;
		//	}
		//}
	}


};

#endif /* GATB_SIMKA_SRC_SUBSAMPLING_SUBSAMPLING_HPP_ */
