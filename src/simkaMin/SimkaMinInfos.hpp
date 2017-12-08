/*
 * SimkaMinInfos.hpp
 *
 *  Created on: 7 déc. 2017
 *      Author: gbenoit
 */

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMININFOS_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMININFOS_HPP_


#include "SimkaMinCommons.hpp"

class SimkaMinInfosAlgorithm : public Algorithm
{

public:

	IProperties* _options;
	string _inputFilename;
	u_int32_t _nbDatasets;
	u_int32_t _sketchSize;



	SimkaMinInfosAlgorithm(IProperties* options):
		Algorithm("simkaMinInfosAlgorithm", -1, options)
	{
	}

	void execute(){

		parseArgs();

		printInfos();

	}

	void parseArgs(){

		_options = getInput();

		_inputFilename = _options->getStr(STR_URI_INPUT);

		if(!System::file().doesExist(_inputFilename)){
			std::cerr << "Error: input does not exist (" << _inputFilename << ")" << std::endl;
			exit(1);
		}
	}

	void printInfos(){

		//vector<string> datasetIds;
		//SimkaMinCommons::readIds(_inputFilename, datasetIds);

		u_int32_t seed;
		u_int8_t kmerSize;
		SimkaMinCommons::getKmerInfos(_inputFilename, kmerSize, _sketchSize, seed, _nbDatasets);

		cout << "Sketch info: " << _inputFilename << endl;

		cout << endl;
		cout << "k-mer size  : " << (u_int32_t) kmerSize << endl;
		cout << "Sketch size : " << _sketchSize << endl;
		cout << "Seed        : " << seed << endl;

		cout << endl;
		cout << "Nb Datasets: " << _nbDatasets << endl;
		printIds();

		cout << endl;
	}

	void printIds(){

		ifstream file(_inputFilename.c_str(), ios::binary);
		file.seekg(SimkaMinCommons::getFilePosition_sketchIds(_nbDatasets, _sketchSize));

		//u_int32_t nbDatasets;
		//file.read((char*)(&nbDatasets), sizeof(nbDatasets));
		string datasetId;

		for(size_t i=0; i<_nbDatasets; i++){
			SimkaMinCommons::readString(datasetId, file);
			cout << datasetId << endl;
			//datasetIds.push_back(datasetId);
		}

		file.close();

	}



};

















class SimkaMinInfos : public Tool{
public:


	SimkaMinInfos(): Tool ("SimkaMin-Infos"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "filename to a sketch file", true));
		parser->getParser (STR_NB_CORES)->setVisible (false);
		parser->getParser (STR_VERBOSE)->setVisible (false);

	}


	void execute ()
	{
		IProperties* args = getInput();

		SimkaMinInfosAlgorithm* algo = new SimkaMinInfosAlgorithm(args);
		algo->execute();
		delete algo;
	}

};

#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMININFOS_HPP_ */







































