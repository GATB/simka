/*
 * SimkaMinDistance.hpp
 *
 *  Created on: 26 mai 2017
 *      Author: gbenoit
 */

#ifndef SIMKA2_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_
#define SIMKA2_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_


template<class CountType>
class SimkaMinDistance {

	size_t _nbBanks;
    vector<size_t> _sharedBanks;

	vector<vector<u_int64_t> > _brayCurtisNumerator;
	vector<u_int64_t> _nbKmersPerDataset;

	vector<vector<float> > _distanceMatrix;

public:

	SimkaMinDistance(){
	}

	SimkaMinDistance(size_t nbBanks){
		_nbBanks = nbBanks;

		_nbKmersPerDataset.resize(_nbBanks);

		//u_int64_t symetricDistanceMatrixSize = (_nbBanks*(_nbBanks+1))/2;
		//_brayCurtisNumerator.resize(symetricDistanceMatrixSize, 0);
		//vector<vector<u_int64_t> > _matrix_squaredHalf;

		_brayCurtisNumerator.resize(_nbBanks-1);

		u_int64_t size = _nbBanks-1;
		for(size_t i=0; i<_brayCurtisNumerator.size(); i++){
			_brayCurtisNumerator[i].resize(size);
			size -= 1;
		}


	}

	~SimkaMinDistance(){

	}

	void processAbundanceVector(vector<CountType>& counts){
		//void processAbundanceVector(vector<KmerCountType>& counts){
		//	for(size_t i=0; i<counts.size(); i++){
		//		cout << counts[i] << " ";
		//	}
		//	cout << endl;
		//}

		_sharedBanks.clear();

		for(size_t i=0; i<counts.size(); i++){
			if(counts[i] > 0){
				//_stats->_nbDistinctKmersPerDataset[i] += 1;
				_nbKmersPerDataset[i] += counts[i];
				_sharedBanks.push_back(i);
			}
		}

		for(size_t ii=0; ii<_sharedBanks.size(); ii++){

			u_int64_t i = _sharedBanks[ii];
			u_int64_t abundanceI = counts[i];
			u_int64_t jOffset = _brayCurtisNumerator.size() - _brayCurtisNumerator[i].size();// + 1;

			for(size_t jj=ii+1; jj<_sharedBanks.size(); jj++){

				u_int64_t j = _sharedBanks[jj];

				u_int64_t abundanceJ = counts[j];

				_brayCurtisNumerator[i][j-jOffset-1] += min(abundanceI, abundanceJ);

			}
		}

		/*
		for(size_t ii=0; ii<_sharedBanks.size(); ii++){
			for(size_t jj=ii+1; jj<_sharedBanks.size(); jj++){

				u_int16_t i = _sharedBanks[ii];
				u_int16_t j = _sharedBanks[jj];
				size_t symetricIndex = j + ((_nbBanks-1)*i) - (i*(i-1)/2);

				u_int64_t abundanceI = counts[i];
				u_int64_t abundanceJ = counts[j];

				//_stats->_matrixNbSharedKmers[i][j] += counts[i];
				//_stats->_matrixNbSharedKmers[j][i] += counts[j];
				//_stats->_matrixNbDistinctSharedKmers[symetricIndex] += 1;

				//cout << i << " " << j << "    " << (j + ((_nbBanks-1)*i) - (i*(i-1)/2)) << endl;
				_brayCurtisNumerator[symetricIndex] += min(abundanceI, abundanceJ);
			}
		}*/
	}

	void computeDistanceMatrix(const string& outputDir, const string& distanceName){
		_distanceMatrix.resize(_nbBanks-1);

		u_int64_t size = _nbBanks-1;
		for(size_t i=0; i<_distanceMatrix.size(); i++){
			_distanceMatrix[i].resize(size, 0);
			size -= 1;
		}




		for(size_t i=0; i<_brayCurtisNumerator.size(); i++){

			u_int64_t jOffset = _brayCurtisNumerator.size() - _brayCurtisNumerator[i].size();

			for(size_t j=0; j<_brayCurtisNumerator[i].size(); j++){

				double dist = distance_abundance_brayCurtis(i, j+1+jOffset, i, j);
				_distanceMatrix[i][j] = dist;
			}
		}

		outputDistanceMatrix(outputDir, distanceName);
	}



    double distance_abundance_brayCurtis(size_t i, size_t j, size_t i2, size_t j2){

    	//double intersection = _stats._abundance_jaccard_intersection[i][j];
    	double union_ = _nbKmersPerDataset[i] + _nbKmersPerDataset[j];

    	if(union_ == 0) return 1;

    	double intersection = 2 * _brayCurtisNumerator[i2][j2];

    	double jaccard = 1 - intersection / union_;

    	return jaccard;
    }


    void outputDistanceMatrix(const string& outputDir, const string& distanceName){

		string filename = outputDir + "/" + distanceName + ".bin";

		ofstream outputFile(filename.c_str(), ios::binary);


		//u_int64_t nbNewBanks = _brayCurtisNumerator.size() + 1;
		//u_int64_t nbBanks = nbOldBanks + nbNewBanks;

		for(size_t i=0; i<_nbBanks; i++){
			writeMatrixBinaryFromSplits_squaredHalf(outputFile, i);
		}

		outputFile.close();

    }

    void writeMatrixBinaryFromSplits_squaredHalf(ofstream& outputFile, size_t i){
		//u_int64_t nbNewBanks = distanceMatrix_squaredHalf.size() + 1;

		//for(size_t i=0; i<nbNewBanks; i++){

		u_int64_t jOffset = _distanceMatrix.size() - _distanceMatrix[i].size();

		for(size_t j=0; j<_nbBanks; j++){

			//cout << i << " " << j << endl;
			if(i == j){
				float distanceZero = 0;
				//cout << distanceZero << "\t\t";
				outputFile.write((const char*)&distanceZero, sizeof(distanceZero));
			}
			else if(j < i){

				u_int64_t iOffset = _distanceMatrix.size() - _distanceMatrix[j].size();

				float value = _distanceMatrix[j][i-iOffset-1];
				//cout << value << "\t\t";
				outputFile.write((const char*)&value, sizeof(value));
			}
			else{
				float value = _distanceMatrix[i][j-jOffset-1];
				//cout << value << "\t\t";
				outputFile.write((const char*)&value, sizeof(value));
			}
			//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._brayCurtisNumerator._matrix_squaredHalf);
			//_matrix_squaredHalf[i][j] = dist;
		}

	}

    void writeMatrixASCII(const string& outputDir, const string& outputDirTemp, const string& distanceName, const vector<string>& bankNames){

		vector<float> rowData(bankNames.size(), 0);
		string filenameBin = outputDirTemp + "/" + distanceName + ".bin";
		ifstream binaryMatrixFile(filenameBin.c_str(), ios::binary);
		string filename = outputDir + "/" + distanceName + ".csv";
		gzFile out = gzopen((filename + ".gz").c_str(),"wb");

		string str = "";

		for(size_t i=0; i<_nbBanks; i++){
			str += ";" + bankNames[i];
		}
		str += '\n';
		gzwrite(out, str.c_str(), str.size());


		for(size_t i=0; i<_nbBanks; i++){

			str = "";
			str += bankNames[i]+ ";";

			//size_t rowIndex = _wantedIdsIndex[i];
			loadRow(i, binaryMatrixFile, rowData);

			for(size_t j=0; j<_nbBanks; j++){

				str += Stringify::format("%f", rowData[j]) + ";";

			}

			str.erase(str.size()-1);
			str += '\n';

			gzwrite(out, str.c_str(), str.size());
		}


		gzclose(out);
		binaryMatrixFile.close();
	}

	void loadRow(size_t rowIndex, ifstream& matrixBinaryFile, vector<float>& resultRow){

		matrixBinaryFile.seekg(rowIndex*resultRow.size()*sizeof(float), ios_base::beg);

		matrixBinaryFile.read((char*)resultRow.data(), sizeof(float)*resultRow.size());

	}

};



#endif /* SIMKA2_SRC_SIMKAMIN_SIMKAMINDISTANCE_HPP_ */
