/*****************************************************************************
 *   SimkaMin: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2019  INRIA
 *   Authors: G.Benoit
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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_


#include <unordered_map>
#include <gatb/gatb_core.hpp>

#define KMER_SPECTRUM_HEADER_SIZE (1+4+4+4) //At the begining of the .kmers file we store the size of the kmer (on 1 byte), the sketch size (on 4 bytes), the seed used by Murmurhash3 (4 bytes), the number of datasets in the sketch file (4 bytes)

const string STR_SIMKA_SEED = "-seed";
const string STR_SIMKA_SKETCH_SIZE = "-nb-kmers";
const string STR_SIMKA_URI_INPUT_1 = "-in1";
const string STR_SIMKA_URI_INPUT_2 = "-in2";
const string STR_SIMKA_INPUT_IDS = "-in-ids";
const string STR_SIMKA_ABUNDANCE_FILTER = "-filter";
//const string STR_SIMKA2_DATASET_ID = "-id";

typedef u_int32_t KmerCountType;
typedef unordered_map<u_int64_t, KmerCountType> KmerCountDictionaryType;
typedef float DistanceValueType;

struct KmerAndCountType{

public:
	u_int64_t _kmer;
	KmerCountType _count;

	KmerAndCountType(){
	}

	KmerAndCountType(u_int64_t kmer, KmerCountType count){
		_kmer = kmer;
		_count = count;
	}
};

struct PairwiseDistance{
	u_int64_t _i;
	u_int64_t _j;
	DistanceValueType _distance;

	PairwiseDistance(){
		_i = -1;
		_j = -1;
		_distance = -1;
	}

	void set(u_int64_t i, u_int64_t j, DistanceValueType distance){
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

		u_int8_t kmerSize;
		u_int32_t sketchSize, seed, nbDatasets;
		getKmerInfos(filename, kmerSize, sketchSize, seed, nbDatasets);

		ifstream file(filename.c_str(), ios::binary);
		file.seekg(SimkaMinCommons::getFilePosition_sketchIds(nbDatasets, sketchSize));

		//u_int32_t nbDatasets;
		//file.read((char*)(&nbDatasets), sizeof(nbDatasets));
		string datasetId;

		for(size_t i=0; i<nbDatasets; i++){
			SimkaMinCommons::readString(datasetId, file);
			datasetIds.push_back(datasetId);
		}

		file.close();
	}

	/*
	static u_int32_t readNbDatasets(const string& filename){

		string filenameIds = filename + ".ids";
		ifstream file(filenameIds.c_str(), ios::binary);

		u_int32_t nbDatasets;
		file.read((char*)(&nbDatasets), sizeof(nbDatasets));

		file.close();

		return nbDatasets;
	}*/

	static void getKmerInfos(const string& filename, u_int8_t& kmerSize, u_int32_t& sketchSize, u_int32_t& seed, u_int32_t& nbDatasets){

		//string filenameKmers = filename + ".kmers";
		ifstream file(filename.c_str(), ios::binary);

		u_int8_t kmerSize_;
		file.read((char*)(&kmerSize_), sizeof(kmerSize_));
		u_int32_t sketchSize_;
		file.read((char*)(&sketchSize_), sizeof(sketchSize_));
		u_int32_t seed_;
		file.read((char*)(&seed_), sizeof(seed_));
		u_int32_t nbDatasets_;
		file.read((char*)(&nbDatasets_), sizeof(nbDatasets_));

		file.close();

		kmerSize = kmerSize_;
		sketchSize = sketchSize_;
		seed = seed_;
		nbDatasets = nbDatasets_;

	}

	static u_int64_t getFilePosition_sketchIds(u_int32_t nbDatasets, u_int32_t sketchSize){
		u_int64_t filePos = KMER_SPECTRUM_HEADER_SIZE + (nbDatasets * sketchSize * sizeof(KmerAndCountType));
		return filePos;
	}

	static u_int64_t getFilePosition_nbDatasets(){
		return 1+4+4;
	}

};









#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINCOMMONS_HPP_ */
