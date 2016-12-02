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



	static void writeMatrixBinaryFromSplits(const string& distanceMatricesDir, const string& distanceName, const vector<vector<float> >& distanceMatrix_rectangular, const vector<vector<float> >& distanceMatrix_squaredHalf){

		/*

		cout << endl;
		cout << endl;
		for(size_t i=0; i<distanceMatrix_squaredHalf.size(); i++){
			for(size_t j=0; j<distanceMatrix_squaredHalf[i].size(); j++){
				cout << distanceMatrix_squaredHalf[i][j] << " ";
			}
			cout << endl;
		}*/

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

		u_int64_t nbOldBanks = 0;
		if(distanceMatrix_rectangular.size() > 0){
			nbOldBanks = distanceMatrix_rectangular[0].size();
		}


		u_int64_t nbNewBanks = distanceMatrix_squaredHalf.size() + 1;
		u_int64_t nbBanks = nbOldBanks + nbNewBanks;

		if(nbOldBanks > 0){
			if(nbNewBanks > 1){
				for(size_t i=0; i<nbNewBanks; i++){
					writeMatrixBinaryFromSplits_rectangular(outputFile, i, distanceMatrix_rectangular);
					writeMatrixBinaryFromSplits_squaredHalf(outputFile, i, distanceMatrix_squaredHalf);
					//cout << endl;
				}
			}
			else{
				for(size_t i=0; i<nbNewBanks; i++){
					writeMatrixBinaryFromSplits_rectangular(outputFile, i, distanceMatrix_rectangular);
				}
			}
		}
		else{
			for(size_t i=0; i<nbNewBanks; i++){
				writeMatrixBinaryFromSplits_squaredHalf(outputFile, i, distanceMatrix_squaredHalf);
			}
		}


		//cout << distanceMatrix.size() << "   " << distanceMatrix[0].size() << endl;

		//for(size_t i=0; i<distanceMatrix.size(); i++){
		//	outputFile.write((const char*)distanceMatrix[i].data(), sizeof(float)*distanceMatrix[i].size());
		//}

		outputFile.close();
	}

	static void writeMatrixBinaryFromSplits_rectangular(ofstream& outputFile, size_t i, const vector<vector<float> >& distanceMatrix_rectangular){

		//cout << endl;
		//cout << distanceMatrix_rectangular.size() << " " << distanceMatrix_rectangular[i].size() << endl;
		//cout << endl;
		//for(size_t j=0; j<distanceMatrix_rectangular[i].size(); j++){
		//	cout << distanceMatrix_rectangular[i][j] << "\t";
		//}

		outputFile.write((const char*)distanceMatrix_rectangular[i].data(), sizeof(float)*distanceMatrix_rectangular[i].size());
	}

	static void writeMatrixBinaryFromSplits_squaredHalf(ofstream& outputFile, size_t i, const vector<vector<float> >& distanceMatrix_squaredHalf){
		u_int64_t nbNewBanks = distanceMatrix_squaredHalf.size() + 1;

		//for(size_t i=0; i<nbNewBanks; i++){

			u_int64_t jOffset = distanceMatrix_squaredHalf.size() - distanceMatrix_squaredHalf[i].size();

			for(size_t j=0; j<nbNewBanks; j++){

				//cout << i << " " << j << endl;
				if(i == j){
					float distanceZero = 0;
					//cout << distanceZero << "\t\t";
					outputFile.write((const char*)&distanceZero, sizeof(distanceZero));
				}
				else if(j < i){

					u_int64_t iOffset = distanceMatrix_squaredHalf.size() - distanceMatrix_squaredHalf[j].size();

					float value = distanceMatrix_squaredHalf[j][i-iOffset-1];
					//cout << value << "\t\t";
					outputFile.write((const char*)&value, sizeof(value));
				}
				else{
					float value = distanceMatrix_squaredHalf[i][j-jOffset-1];
					//cout << value << "\t\t";
					outputFile.write((const char*)&value, sizeof(value));
				}
				//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._brayCurtisNumerator._matrix_squaredHalf);
				//_matrix_squaredHalf[i][j] = dist;
			}

			//cout << endl;
		//}
	}

	static void writeMatrixBinary(const string& distanceMatricesDir, const string& distanceName, const vector<vector<float> >& distanceMatrix){
		string filename = distanceMatricesDir + "/" + distanceName + ".bin";

		ofstream outputFile(filename.c_str(), ios::binary);

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
