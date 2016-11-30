/*
 * SimkaDistanceMatrixBinary.hpp
 *
 *  Created on: 15 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_CORE_SIMKADISTANCEMATRIXBINARY_HPP_
#define GATB_SIMKA_SRC_CORE_SIMKADISTANCEMATRIXBINARY_HPP_

#include "../utils/SimkaIoUtils.hpp"

class SimkaDistanceMatrixBinary {
public:


	static void saveMatrixIds(const string& dir, const vector<string>& bankNames, u_int64_t nbProcessedDatasets, u_int64_t nbNewDatasets){

		string matrixInfoFilename = dir + "/matrix_infos.bin";
		ofstream matrixInfoFile(matrixInfoFilename.c_str(), std::ios::binary);

		matrixInfoFile.write((char const*)(&nbNewDatasets), sizeof(nbNewDatasets));

		for(size_t i=nbProcessedDatasets; i<(nbProcessedDatasets+nbNewDatasets); i++){
			SimkaIoUtils::SimkaIoUtils::simka2_writeString(bankNames[i], matrixInfoFile);
		}

		matrixInfoFile.close();

	}

	static bool loadMatrixIds(const string& dir, vector<string>& bankNames){
		u_int64_t nbBanks;

		string matrixInfoFilename = dir + "/matrix_infos.bin";
		if(!System::file().doesExist(matrixInfoFilename)){
			bankNames.clear();
			return false;
		}

		ifstream matrixInfoFile(matrixInfoFilename.c_str(), std::ios::binary);

		matrixInfoFile.read((char*)(&nbBanks), sizeof(nbBanks));
		//cout << "nb banks: " << nbBanks << endl;
		bankNames.resize(nbBanks);

		for(size_t i=0; i<bankNames.size(); i++){
			string bankName;
			SimkaIoUtils::simka2_readString(bankName, matrixInfoFile);
			bankNames[i] = bankName;
		}

		matrixInfoFile.close();

		return true;
	}


	static void writeMatrixBinary(const string& distanceMatricesDir, const string& distanceName, const vector<vector<float> >& distanceMatrix){

		/*
		cout << "lala write" << endl;
		for(size_t i=0; i<distanceMatrix.size(); i++){

			for(size_t j=0; j<distanceMatrix[i].size(); j++){
				cout << distanceMatrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;*/
		//string distanceMatrixDir = distanceMatricesDir + "/" + distanceName;
		//if(System::file().doesExist(distanceMatrixDir)){

		//}
		//else{
		//	System::file().mkdir(distanceMatrixDir, -1);
		//}

		//string distanceMatrixDir = outputDirTemp + "/distance_matrix";
		string filename = distanceMatricesDir + "/" + distanceName + ".bin";

		ofstream outputFile(filename.c_str(), ios::binary);

		//cout << distanceMatrix.size() << "   " << distanceMatrix[0].size() << endl;

		for(size_t i=0; i<distanceMatrix.size(); i++){
			outputFile.write((const char*)distanceMatrix[i].data(), sizeof(float)*distanceMatrix[i].size());
		}

		outputFile.close();
	}

	static void loadMatrix(const string& filename, vector<vector<float> >& distanceMatrix){

		//size_t nbBanks = distanceMatrix.size();
		//size_t nbNewBanks = distanceMatrix[0].size();
		//size_t nbOldBanks = distanceMatrix.size();
		//size_t nbBanks = nbOldBanks + nbNewBanks;

		ifstream inputFile(filename.c_str(), ios::binary);

		for(size_t i=0; i<distanceMatrix.size(); i++){
			inputFile.read((char*)distanceMatrix[i].data(), sizeof(float)*distanceMatrix[i].size());
		}

		inputFile.close();

		/*
		cout << "lala read" << endl;
		for(size_t i=0; i<distanceMatrix.size(); i++){

			for(size_t j=0; j<distanceMatrix[i].size(); j++){
				cout << distanceMatrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;*/

	}

	//Result must be of size nbBanks
	static void loadRow(size_t rowIndex, ifstream& matrixBinaryFile, vector<float>& resultRow){

		matrixBinaryFile.seekg(rowIndex*resultRow.size()*sizeof(float), ios_base::beg);

		matrixBinaryFile.read((char*)resultRow.data(), sizeof(float)*resultRow.size());

	}

};

#endif /* GATB_SIMKA_SRC_CORE_SIMKADISTANCEMATRIXBINARY_HPP_ */
