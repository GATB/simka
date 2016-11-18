/*
 * SimkaIoUtils.h
 *
 *  Created on: 18 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_
#define GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_

class SimkaIoUtils {
public:
	//SimkaIoUtils();
	//virtual ~SimkaIoUtils();

	static string getDatasetID(const string& kmerSpectrumDir){
		string datasetID = System::file().getBaseName(kmerSpectrumDir);
		//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
		return datasetID;
	}



	static u_int64_t simka2_getFileSize(const string& filename){
		std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
		u_int64_t size = in.tellg();
		in.close();
		return size;
	}
	/*
	u_int64_t simka2_getDatasetSize(const string& kmerSpectrumDir){
		u_int64_t datasetSize = 0;
		for(size_t partitionID=0; partitionID<SIMKA2_NB_PARTITIONS; partitionID++){
			string partFilename = kmerSpectrumDir + "/" + Stringify::format("%i", partitionID) + ".gz";
			datasetSize += simka2_getFileSize(partFilename);
		}
		return datasetSize;
	}
	*/

	static void simka2_writeString(const string& s, ofstream& file){
		u_int16_t size = s.size();
		file.write((char const*)(&size), sizeof(size));
		file.write(s.c_str(), size);
	}

	static void simka2_readString(string& s, ifstream& file){
		u_int16_t size;
		file.read((char*)(&size), sizeof(size));
		std::vector<char> buffer(size);
		file.read(&buffer[0], buffer.size());
		s.assign(buffer.begin(), buffer.end());
		//return string linkedDatasetID( buffer.begin(), buffer.end() );
	}

	static void simka2_writeDatasetInfo(ofstream& file, const string& datasetID, u_int64_t nbReads, u_int64_t nbDistinctKmers, u_int64_t nbKmers, u_int64_t chord_N2){
	    simka2_writeString(datasetID, file);
	    file.write((char const*)(&nbReads), sizeof(nbReads));
	    file.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	    file.write((char const*)(&nbKmers), sizeof(nbKmers));
	    file.write((char const*)(&chord_N2), sizeof(chord_N2));
	}

	static void simka2_readDatasetInfo(ifstream& file, string& datasetID, u_int64_t& nbReads, u_int64_t& nbDistinctKmers, u_int64_t& nbKmers, u_int64_t& chord_N2){
		simka2_readString(datasetID, file);
		file.read((char *)(&nbReads), sizeof(nbReads));
		file.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
		file.read((char *)(&nbKmers), sizeof(nbKmers));
		file.read((char *)(&chord_N2), sizeof(chord_N2));
	}

	static void simka2_transferDatasetInfo(ifstream& source, ofstream& dest){
		string datasetID;
		u_int64_t nbReads;
		u_int64_t nbDistinctKmers;
		u_int64_t nbKmers;
		u_int64_t chord_N2;

		simka2_readString(datasetID, source);
		source.read((char *)(&nbReads), sizeof(nbReads));
		source.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
		source.read((char *)(&nbKmers), sizeof(nbKmers));
		source.read((char *)(&chord_N2), sizeof(chord_N2));

	    simka2_writeString(datasetID, dest);
	    dest.write((char const*)(&nbReads), sizeof(nbReads));
	    dest.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	    dest.write((char const*)(&nbKmers), sizeof(nbKmers));
	    dest.write((char const*)(&chord_N2), sizeof(chord_N2));
	}

};

#endif /* GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_ */
