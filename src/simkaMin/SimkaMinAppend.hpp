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

#ifndef SIMKA1_4_SRC_SIMKAMIN_SIMKAMINAPPEND_HPP_
#define SIMKA1_4_SRC_SIMKAMIN_SIMKAMINAPPEND_HPP_

#include "SimkaMinCommons.hpp"



/*
Header and sketches of first sketch (-in1) are kept
First sketch overwrite starts after its sketches, where ids starts
Sketches of second files (-in2) are written after sketches of first file (-in1)
Then, Ids of first file are written, and then ids of second file
The number of datasets in the header is updated
 */

class SimkaMinAppendAlgorithm : public Algorithm
{

public:

	IProperties* _options;
	string _inputFilename1;
	string _inputFilename2;
	u_int32_t _nbDatasets;
	u_int32_t _sketchSize;

	ofstream _outputFile;
	ifstream _inputFile2;


	SimkaMinAppendAlgorithm(IProperties* options):
		Algorithm("simkaMinAppendAlgorithm", -1, options)
	{
	}

	void execute(){

		parseArgs();
		append();
	}

	void parseArgs(){

		_options = getInput();

		_inputFilename1 = _options->getStr(STR_SIMKA_URI_INPUT_1);
		_inputFilename2 = _options->getStr(STR_SIMKA_URI_INPUT_2);

		if(!System::file().doesExist(_inputFilename1)){
			std::cerr << "Error: input does not exist (" << _inputFilename1 << ")" << std::endl;
			exit(1);
		}
		if(!System::file().doesExist(_inputFilename2)){
			std::cerr << "Error: input does not exist (" << _inputFilename2 << ")" << std::endl;
			exit(1);
		}
	}

	void append(){

		u_int8_t kmerSize1, kmerSize2;
		u_int32_t sketchSize1, sketchSize2, seed1, seed2, nbDatasets1, nbDatasets2;

		SimkaMinCommons::getKmerInfos(_inputFilename1, kmerSize1, sketchSize1, seed1, nbDatasets1);
		SimkaMinCommons::getKmerInfos(_inputFilename2, kmerSize2, sketchSize2, seed2, nbDatasets2);

		if(kmerSize1 != kmerSize2){
			std::cerr << "Error: can't merge sketches with different kmer sizes (" << kmerSize1 << " vs " << kmerSize2 << ")" << std::endl;
			exit(1);
		}
		if(sketchSize1 != sketchSize2){
			std::cerr << "Error: can't merge sketches with different sketch sizes (" << sketchSize1 << " vs " << sketchSize2 << ")" << std::endl;
			exit(1);
		}
		if(seed1 != seed2){
			std::cerr << "Error: can't merge sketches with different seeds (" << seed1 << " vs " << seed2 << ")" << std::endl;
			exit(1);
		}




		u_int32_t nbDatasets = nbDatasets1 + nbDatasets2;
		vector<string> id1, id2;
		SimkaMinCommons::readIds(_inputFilename1, id1);
		SimkaMinCommons::readIds(_inputFilename2, id2);

		//open first file to be overwritten (but without rewriting all its sketches)
		_outputFile.open(_inputFilename1, ios::binary|ios::in);
		_inputFile2.open(_inputFilename2, ios::binary);

		//Update number of datasets in the header
		_outputFile.seekp(SimkaMinCommons::getFilePosition_nbDatasets());
		_outputFile.write((const char*)&nbDatasets, sizeof(nbDatasets));

		appendSkecthes(nbDatasets1, sketchSize1, nbDatasets2);

		appendIds(id1);
		appendIds(id2);

		_inputFile2.close();
		_outputFile.close();
	}

	void appendSkecthes(u_int32_t nbDatasets1, u_int32_t sketchSize1, u_int32_t nbDatasets2){
		_outputFile.seekp(SimkaMinCommons::getFilePosition_sketchIds(nbDatasets1, sketchSize1));
		_inputFile2.seekg(KMER_SPECTRUM_HEADER_SIZE);


		u_int64_t dataToTransfer = nbDatasets2*sketchSize1*sizeof(KmerAndCountType);

		u_int64_t bufferSize = 1024;
		char  buffer[bufferSize];

		/* copy from input to output */
		while (dataToTransfer > 0) {
			u_int64_t size = min(bufferSize, dataToTransfer);
			_inputFile2.read(buffer, size);
			_outputFile.write(buffer, size);
			dataToTransfer -= size;
		}

		//fclose(infile);
		//fclose(outfile);
	}

	void appendIds(vector<string>& ids){

		for(size_t i=0; i<ids.size(); i++){

			string& bankId = ids[i];

			u_int8_t idSize = bankId.size();
			_outputFile.write((const char*)& idSize, sizeof(idSize));
			_outputFile.write(bankId.c_str(), bankId.size());

		}

	}

};

















class SimkaMinAppend : public Tool{
public:


	SimkaMinAppend(): Tool ("SimkaMin-Append"){

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_2, "second sketch file to merge (this file will be appended to the first one)", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA_URI_INPUT_1, "first sketch file to merge (this file will be overwritten)", true));
		parser->getParser (STR_NB_CORES)->setVisible (false);
		parser->getParser (STR_VERBOSE)->setVisible (false);

	}


	void execute ()
	{
		IProperties* args = getInput();

		SimkaMinAppendAlgorithm* algo = new SimkaMinAppendAlgorithm(args);
		algo->execute();
		delete algo;
	}

};

#endif /* SIMKA1_4_SRC_SIMKAMIN_SIMKAMINAPPEND_HPP_ */
