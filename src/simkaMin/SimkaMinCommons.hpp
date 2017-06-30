/*
 * SimkaMinCommons.hpp
 *
 *  Created on: 24 juin 2017
 *      Author: gbenoit
 */

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_


#include <unordered_map>
#include <gatb/gatb_core.hpp>

#define KMER_SPECTRUM_HEADER_SIZE (1+4) //At the begining of the .kmers file we store the size of the kmer (on 1 byte) and the sketch size (on 4 bytes)

const string STR_SIMKA_SKETCH_SIZE = "-nb-kmers";
const string STR_SIMKA_URI_INPUT_1 = "-in1";
const string STR_SIMKA_URI_INPUT_2 = "-in2";
const string STR_SIMKA_INPUT_IDS = "-in-ids";
//const string STR_SIMKA2_DATASET_ID = "-id";

typedef u_int16_t KmerCountType;
typedef unordered_map<u_int64_t, KmerCountType> KmerCountDictionaryType;
typedef float DistanceValueType;

struct KmerAndCountType{

public:
	u_int64_t _kmer;
	KmerCountType _count;

	KmerAndCountType(u_int64_t kmer, KmerCountType count){
		_kmer = kmer;
		_count = count;
	}
};

struct PairwiseDistance{
	u_int64_t _i;
	u_int64_t _j;
	DistanceValueType _distance;

	PairwiseDistance(u_int64_t i, u_int64_t j, DistanceValueType distance){
		_i = i;
		_j = j;
		_distance = distance;
	}
};

class SimkaMinCommons {
public:
	//SimkaMinCommons();
	//virtual ~SimkaMinCommons();

	static void writeString(const string& s, ofstream& file){
		u_int8_t size = s.size();
		file.write((char const*)(&size), sizeof(size));
		file.write(s.c_str(), size);
	}

	static void readString(string& s, ifstream& file){
		u_int8_t size;
		file.read((char*)(&size), sizeof(size));
		std::vector<char> buffer(size);
		file.read(&buffer[0], buffer.size());
		s.assign(buffer.begin(), buffer.end());
		//return string linkedDatasetID( buffer.begin(), buffer.end() );
	}

	static void readIds(const string& filename, vector<string>& datasetIds){

		string filenameIds = filename + ".ids";
		ifstream file(filenameIds.c_str(), ios::binary);

		u_int32_t nbDatasets;
		file.read((char*)(&nbDatasets), sizeof(nbDatasets));
		string datasetId;

		for(size_t i=0; i<nbDatasets; i++){
			readString(datasetId, file);
			datasetIds.push_back(datasetId);
		}

		file.close();
	}

	static u_int32_t readNbDatasets(const string& filename){

		string filenameIds = filename + ".ids";
		ifstream file(filenameIds.c_str(), ios::binary);

		u_int32_t nbDatasets;
		file.read((char*)(&nbDatasets), sizeof(nbDatasets));

		file.close();

		return nbDatasets;
	}

	static void getKmerInfos(const string& filename, size_t& kmerSize, size_t& sketchSize){

		string filenameKmers = filename + ".kmers";
		ifstream file(filenameKmers.c_str(), ios::binary);

		u_int8_t kmerSize_;
		file.read((char*)(&kmerSize_), sizeof(kmerSize_));
		u_int32_t sketchSize_;
		file.read((char*)(&sketchSize_), sizeof(sketchSize_));

		file.close();

		kmerSize = kmerSize_;
		sketchSize = sketchSize_;
	}


};









#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_ */
