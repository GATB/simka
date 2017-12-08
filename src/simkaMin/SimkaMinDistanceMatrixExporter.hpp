/*
 * SimkaMinDistanceMatrixExporter.h
 *
 *  Created on: 29 juin 2017
 *      Author: gbenoit
 */

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXEXPORTER_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXEXPORTER_HPP_



#include "SimkaMinCommons.hpp"






class SimkaDistanceMatrixBinary {
public:

	static void loadRow(size_t rowIndex, ifstream& matrixBinaryFile, vector<float>& resultRow){

		matrixBinaryFile.seekg(rowIndex*resultRow.size()*sizeof(float), ios_base::beg);

		matrixBinaryFile.read((char*)resultRow.data(), sizeof(float)*resultRow.size());

	}

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
			string matrixFilename = _inputDir + "/" + matrixFilenames[i];
			if(matrixFilename.find("mat_") == string::npos) continue;
			if(matrixFilename.find(".bin") == string::npos) continue;

			string distanceName = matrixFilenames[i];
			string ext = ".bin";
			std::string::size_type it = distanceName.find(ext);
			distanceName.erase(it, ext.length());

			cout << distanceName << endl;
			writeMatrixASCII(distanceName, matrixFilename);
		}

	}

	void writeMatrixASCII(const string& distanceName, const string& binaryMatrixFilename){

		vector<float> rowData(_ids1.size(), 0);
		ifstream binaryMatrixFile(binaryMatrixFilename.c_str(), ios::binary);
		string filename = _outputDir + "/" + distanceName + ".csv";
		gzFile out = gzopen((filename + ".gz").c_str(),"wb");

		string str = "";

		for(size_t i=0; i<_ids1.size(); i++){
			str += ";" + _ids1[i]; //_ids[_wantedIdsIndex[i]];
		}
		str += '\n';
		gzwrite(out, str.c_str(), str.size());

		for(size_t i=0; i<_ids2.size(); i++){

			str = "";
			str += _ids2[i] + ";"; //[_wantedIdsIndex[i]] + ";";

			//size_t rowIndex = _wantedIdsIndex[i];
			SimkaDistanceMatrixBinary::loadRow(i, binaryMatrixFile, rowData);

			for(size_t j=0; j<_ids1.size(); j++){

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
