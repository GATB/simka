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







































