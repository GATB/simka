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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXMERGER_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXMERGER_HPP_



#include "SimkaMinCommons.hpp"
#include "SimkaMinDistanceMatrixExporter.hpp"


class SimkaMinDistanceMatrixMergerAlgorithm : public Algorithm{
public:


	IProperties* _options;
	string _inputDir;
	//string _outputDir;
	string _inputSketchFilename_existingDatasets;
	string _inputSketchFilename_newDatasets;

	//vector<string> _ids1;
	//vector<string> _ids2;
	//vector<string> _wantedIds;
	//vector<size_t> _wantedIdsIndex_1;
	//vector<size_t> _wantedIdsIndex_2;
	//unordered_map<string, size_t> _idToIndex_1;
	//unordered_map<string, size_t> _idToIndex_2;
	//size_t _inputMatrixSize_1;
	//size_t _inputMatrixSize_2;
	//size_t _outputMatrixSize;


	SimkaMinDistanceMatrixMergerAlgorithm(IProperties* options):
		Algorithm("simkaMinDistanceMatrixMergerAlgorithm", -1, options)
	{
	}

	void execute(){

		parseArgs();
		mergeMatrices();
	}

	void parseArgs(){

		_options = getInput();

		_inputDir = _options->getStr(STR_URI_INPUT) + "/";
		_inputSketchFilename_existingDatasets = _options->getStr(STR_SIMKA_URI_INPUT_1);
		_inputSketchFilename_newDatasets = _options->getStr(STR_SIMKA_URI_INPUT_2);
		//_outputDir = _options->getStr(STR_URI_OUTPUT);

		//if(getInput()->get(STR_SIMKA_INPUT_IDS)){
		//	_inputFilenameIds =  getInput()->getStr(STR_SIMKA_INPUT_IDS);
		//}

		//if(!System::file().doesExist(_outputDir)){
		//	int ok = System::file().mkdir(_outputDir, -1);
		//	if(ok != 0){
		//		  std::cerr << "Error: can't create output directory (" << _outputDir << ")" << std::endl;
		//		  exit(1);
		//	}
		//}
	}




	void mergeMatrices(){

		u_int32_t dummy;
		u_int8_t dummy_k;
		u_int32_t nbDatasets_existing, nbDatasets_new;
		SimkaMinCommons::getKmerInfos(_inputSketchFilename_existingDatasets, dummy_k, dummy, dummy, nbDatasets_existing);
		SimkaMinCommons::getKmerInfos(_inputSketchFilename_newDatasets, dummy_k, dummy, dummy, nbDatasets_new);


		vector<string> matrixFilenames = System::file().listdir(_inputDir);

		for(size_t i=0; i<matrixFilenames.size(); i++){
			string matrixFilename = matrixFilenames[i];
			if(matrixFilename.find("mat_") == string::npos) continue;
			if(matrixFilename.find(".bin") == string::npos) continue;
			if(matrixFilename.find(".temp") != string::npos) continue;
			//matrixFilename = _inputDir + "/" + matrixFilename;

			cout << matrixFilename << endl;

			//string ext = ".bin";
			//std::string::size_type it = distanceName.find(ext);
			//distanceName.erase(it, ext.length());

			//static void mergeMatrices(const string& existingMatrixFilename, const string& newMatrixFilename_existingVsNew, const string& newMatrixFilename_newVsNew, u_int32_t nbDatasets_existing, u_int32_t nbDatasets_new){

			//cout << distanceName << endl;
			string matrixFilename_existing = _inputDir + matrixFilename;
			string newMatrixFilename_existingVsNew = _inputDir + "/existingVsNew/" + matrixFilename;
			string newMatrixFilename_newVsNew = _inputDir + "/newVsNew/" + matrixFilename;
			SimkaDistanceMatrixBinary::mergeMatrices(matrixFilename_existing, newMatrixFilename_existingVsNew, newMatrixFilename_newVsNew, nbDatasets_existing, nbDatasets_new);
			System::file().remove(matrixFilename_existing);
			System::file().rename(matrixFilename_existing + ".temp", matrixFilename_existing);
			System::file().remove(matrixFilename_existing + ".temp");
		}


	}




};







class SimkaMinDistanceMatrixMerger : public Tool{
public:


	SimkaMinDistanceMatrixMerger(): Tool ("SimkaMin-DistanceMatrixMerger"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for distance matrices", false, "./simkaMin_results"));
	    //parser->push_front (new OptionOneParam (STR_SIMKA_INPUT_IDS, "filename of ids in the result matrix (one id per line). Do not used this option to used all ids.", false));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_2, "sketch file of new datasets", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_1, "sketch file of existing datasets", true));
	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input dir containing existing simka results", true));

	}


	void execute (){

		IProperties* args = getInput();

		SimkaMinDistanceMatrixMergerAlgorithm* algo = new SimkaMinDistanceMatrixMergerAlgorithm(args);
		algo->execute();
		delete algo;
	}

};


#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINDISTANCEMATRIXMERGER_HPP_ */
