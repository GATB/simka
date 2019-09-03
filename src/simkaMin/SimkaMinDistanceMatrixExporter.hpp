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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXEXPORTER_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXEXPORTER_HPP_



#include "SimkaMinCommons.hpp"






class SimkaDistanceMatrixBinary {
public:

	static void loadRow(size_t rowIndex, ifstream& matrixBinaryFile, vector<float>& resultRow){

		matrixBinaryFile.seekg(rowIndex*resultRow.size()*sizeof(float), ios_base::beg);

		matrixBinaryFile.read((char*)resultRow.data(), sizeof(float)*resultRow.size());

	}

	static void mergeMatrices(const string& existingMatrixFilename, const string& newMatrixFilename_existingVsNew, const string& newMatrixFilename_newVsNew, u_int32_t nbDatasets_existing, u_int32_t nbDatasets_new){

		ifstream existingMatrixFile;
		existingMatrixFile.open(existingMatrixFilename.c_str(), ios::binary);
		vector<float> existingRowData(nbDatasets_existing, 0);

		ifstream matrixFile_existingVsNew;
		matrixFile_existingVsNew.open(newMatrixFilename_existingVsNew.c_str(), ios::binary);
		ifstream matrixFile_newVsNew;
		matrixFile_newVsNew.open(newMatrixFilename_newVsNew.c_str(), ios::binary);
		vector<float> newRowData(nbDatasets_new, 0);

		string tempOutputFilename = existingMatrixFilename + ".temp";
		ofstream tempOutputFile;
		tempOutputFile.open(tempOutputFilename.c_str(), ios::binary);

		//Write existing distance + matrixFile_existingVsNew (right part)
		for(size_t i=0; i<nbDatasets_existing; i++){

			SimkaDistanceMatrixBinary::loadRow(i, existingMatrixFile, existingRowData);
			for(size_t j=0; j<existingRowData.size(); j++){
				float distanceValue = existingRowData[j];
				tempOutputFile.write((const char*)&distanceValue, sizeof(distanceValue));
			}

			SimkaDistanceMatrixBinary::loadRow(i, matrixFile_existingVsNew, newRowData);
			for(size_t j=0; j<newRowData.size(); j++){
				float distanceValue = newRowData[j];
				tempOutputFile.write((const char*)&distanceValue, sizeof(distanceValue));
			}

		}

		matrixFile_existingVsNew.seekg(0);
		u_int64_t rowSize = (nbDatasets_existing+nbDatasets_new) * sizeof(float);
		u_int64_t startPos = nbDatasets_existing * rowSize;


		//Write matrixFile_existingVsNew (bottom part)
		for(size_t i=0; i<nbDatasets_existing; i++){

			//SimkaDistanceMatrixBinary::loadRow(i, existingMatrixFile, existingRowData);
			//for(size_t j=0; j<existingRowData.size(); j++){
			//	float distanceValue = existingRowData[j];
			//	tempOutputFile.write((const char*)&distanceValue, sizeof(distanceValue));
			//}

			SimkaDistanceMatrixBinary::loadRow(i, matrixFile_existingVsNew, newRowData);
			for(size_t j=0; j<newRowData.size(); j++){
				float distanceValue = newRowData[j];
				u_int64_t pos = startPos + (j*rowSize) + (i*sizeof(float));
				tempOutputFile.seekp(pos);
				tempOutputFile.write((const char*)&distanceValue, sizeof(distanceValue));
			}

		}


		//Write matrixFile_newVsNew
		for(size_t i=0; i<nbDatasets_new; i++){

			SimkaDistanceMatrixBinary::loadRow(i, matrixFile_newVsNew, newRowData);
			u_int64_t pos = startPos + (nbDatasets_existing*sizeof(float)) + (i*rowSize);
			tempOutputFile.seekp(pos);

			for(size_t j=0; j<newRowData.size(); j++){
				float distanceValue = newRowData[j];
				tempOutputFile.write((const char*)&distanceValue, sizeof(distanceValue));
			}

		}


		existingMatrixFile.close();
		matrixFile_existingVsNew.close();
		matrixFile_newVsNew.close();
		tempOutputFile.close();
	}

/*
	static void writeMatrixBinaryFromSplits(const string& distanceMatricesDir, const string& distanceName, const vector<vector<float> >& distanceMatrix_rectangular, const vector<vector<float> >& distanceMatrix_squaredHalf){


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
}*/

};





class SimkaMinDistanceMatrixExporterAlgorithm : public Algorithm{
public:


	IProperties* _options;
	string _inputDir;
	string _outputDir;
	string _inputFilenameIds;
	string _inputSketchFilename1;
	string _inputSketchFilename2;

	vector<string> _ids1;
	vector<string> _ids2;
	//vector<string> _wantedIds;
	//vector<size_t> _wantedIdsIndex_1;
	//vector<size_t> _wantedIdsIndex_2;
	//unordered_map<string, size_t> _idToIndex_1;
	//unordered_map<string, size_t> _idToIndex_2;
	size_t _inputMatrixSize_1;
	size_t _inputMatrixSize_2;
	//size_t _outputMatrixSize;


	SimkaMinDistanceMatrixExporterAlgorithm(IProperties* options):
		Algorithm("simkaMinDistanceExporterAlgorithm", -1, options)
	{
	}

	void execute(){
		_inputFilenameIds = "";
		parseArgs();
		createWantedIds();
		//createIdsIndex();
		writeMatrices();
	}

	void parseArgs(){

		_options = getInput();

		_inputDir = _options->getStr(STR_URI_INPUT);
		_inputSketchFilename1 = _options->getStr(STR_SIMKA_URI_INPUT_1);
		_inputSketchFilename2 = _options->getStr(STR_SIMKA_URI_INPUT_2);
		_outputDir = _options->getStr(STR_URI_OUTPUT);

		if(getInput()->get(STR_SIMKA_INPUT_IDS)){
			_inputFilenameIds =  getInput()->getStr(STR_SIMKA_INPUT_IDS);
		}

		if(!System::file().doesExist(_outputDir)){
			int ok = System::file().mkdir(_outputDir, -1);
			if(ok != 0){
				  std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
				  exit(1);
			}
		}
	}


	void createWantedIds(){

		SimkaMinCommons::readIds(_inputSketchFilename1, _ids1);
		SimkaMinCommons::readIds(_inputSketchFilename2, _ids2);

		/*
		if(_inputFilenameIds.empty()){
			_wantedIds = vector<string>(_ids1);
			_wantedIds.insert(_wantedIds.end(), _ids2.begin(), _ids2.end());
		}
		else{
			string line;
			ifstream inputFile(_inputFilenameIds.c_str());

			while(getline(inputFile, line)){

				line.erase(std::remove(line.begin(),line.end(),' '),line.end());
				if(line == "") continue;

				_wantedIds.push_back(line);
			}
		}
		*/

		_inputMatrixSize_1 = _ids1.size();
		_inputMatrixSize_2 = _ids2.size();
		//_outputMatrixSize = _wantedIds.size();
		cout << "Matrix size: " << _inputMatrixSize_1 << " x " << _inputMatrixSize_2 << endl;
	}

	/*
	void createIdsIndex(){

		for(size_t i=0; i<_ids1.size(); i++){
			_idToIndex_1[_ids1[i]] = i;
		}
		for(size_t i=0; i<_ids2.size(); i++){
			_idToIndex_2[_ids2[i]] = i;
		}

		for(size_t i=0; i<_wantedIds.size(); i++){
			if(_idToIndex_1.find(_wantedIds[i]) == _idToIndex_1.end()){
				cout << "ID not found in distance matrix: " << _wantedIds[i] << endl;
			}
			else{
				_wantedIdsIndex.push_back(_idToIndex[_wantedIds[i]]);
			}
		}

		//_wantedIdsIndex.resize(_outputMatrixSize);
		//for(size_t i=0; i<_outputMatrixSize; i++){
		//}
		_outputMatrixSize = _wantedIdsIndex.size();
		cout << "output matrix size: " << _outputMatrixSize << endl;
	}
	*/

	void writeMatrices(){
		vector<string> matrixFilenames = System::file().listdir(_inputDir);


		for(size_t i=0; i<matrixFilenames.size(); i++){
			string matrixFilename = matrixFilenames[i];
			if(matrixFilename.find("mat_") == string::npos) continue;
			if(matrixFilename.find(".bin") == string::npos) continue;
			if(matrixFilename.find(".temp") != string::npos) continue;
			matrixFilename = _inputDir + "/" + matrixFilename;

			string distanceName = matrixFilenames[i];
			string ext = ".bin";
			std::string::size_type it = distanceName.find(ext);
			distanceName.erase(it, ext.length());

			cout << distanceName << endl;
			writeMatrixASCII(distanceName, matrixFilename);
		}

	}

	void writeMatrixASCII(const string& distanceName, const string& binaryMatrixFilename){

		vector<float> rowData(_ids2.size(), 0);
		ifstream binaryMatrixFile(binaryMatrixFilename.c_str(), ios::binary);
		string filename = _outputDir + "/" + distanceName + ".csv";
		gzFile out = gzopen((filename + ".gz").c_str(),"wb");

		string str = "";

		for(size_t i=0; i<_ids2.size(); i++){
			str += ";" + _ids2[i]; //_ids[_wantedIdsIndex[i]];
		}
		str += '\n';
		gzwrite(out, str.c_str(), str.size());

		for(size_t i=0; i<_ids1.size(); i++){

			str = "";
			str += _ids1[i] + ";"; //[_wantedIdsIndex[i]] + ";";

			//size_t rowIndex = _wantedIdsIndex[i];
			SimkaDistanceMatrixBinary::loadRow(i, binaryMatrixFile, rowData);

			for(size_t j=0; j<_ids2.size(); j++){

				//str += Stringify::format("%f", rowData[_wantedIdsIndex[j]]) + ";";
				str += Stringify::format("%f", rowData[j]) + ";";

			}

			str.erase(str.size()-1);
			str += '\n';

			gzwrite(out, str.c_str(), str.size());
		}

		gzclose(out);
		binaryMatrixFile.close();
	}


};







class SimkaMinDistanceMatrixExporter : public Tool{
public:


	SimkaMinDistanceMatrixExporter(): Tool ("SimkaMin-DistanceMatrixExporter"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for distance matrices", false, "./simkaMin_results"));
	    //parser->push_front (new OptionOneParam (STR_SIMKA_INPUT_IDS, "filename of ids in the result matrix (one id per line). Do not used this option to used all ids.", false));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_2, "second used sketch file (-in2 argument of ./simkaMin distance)", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_1, "first used sketch file (-in1 argument of ./simkaMin distance)", true));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input dir containing distance matrices in binary format (-out argument of ./simkaMin distance)", true));

	}


	void execute (){

		IProperties* args = getInput();

		SimkaMinDistanceMatrixExporterAlgorithm* algo = new SimkaMinDistanceMatrixExporterAlgorithm(args);
		algo->execute();
		delete algo;
	}

};


#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXEXPORTER_HPP_ */
