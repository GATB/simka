/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
//#include "../core/SimkaUtils.hpp"
#include "../simka2/Simka2Utils.hpp"
#include "SimkaDistanceMatrixBinary.hpp"
//#include "../core/SimkaUtils.hpp"
//#include "Simka2Database.hpp"
//#include "../minikc/MiniKC.hpp"
//#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
//#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"




/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
class SimkaDistanceExporterAlgorithm : public Algorithm
{

public:

	string _dirBinaryMatrices;
	string _inputFilenameIds;
	string _outputDir;

	vector<string> _ids;
	vector<string> _wantedIds;
	vector<size_t> _wantedIdsIndex;
	map<string, size_t> _idToIndex;
	size_t _matrixSize;

	SimkaDistanceExporterAlgorithm(IProperties* options):
		Algorithm("simkaDistanceExporter", -1, options)
	{
		_inputFilenameIds = "";
	}

	void execute(){
		parseArgs();
		createWantedIds();
		createIdsIndex();
		writeMatrices();
	}

	void parseArgs(){

		_dirBinaryMatrices =  getInput()->getStr(STR_URI_INPUT);
		_outputDir =  getInput()->getStr(STR_URI_OUTPUT);

		if(getInput()->get(STR_SIMKA2_INPUT_IDS)){
			_inputFilenameIds =  getInput()->getStr(STR_SIMKA2_INPUT_IDS);
		}

		if(!System::file().doesExist(_outputDir)){
			System::file().mkdir(_outputDir, -1);
		}

	}

	void createWantedIds(){

		SimkaDistanceMatrixBinary::loadMatrixIds(_dirBinaryMatrices, _ids);

		if(_inputFilenameIds.empty()){
			_wantedIds = vector<string>(_ids);
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

		_matrixSize = _wantedIds.size();
		cout << "matrix size: " << _matrixSize << endl;
	}

	void createIdsIndex(){

		for(size_t i=0; i<_ids.size(); i++){
			_idToIndex[_ids[i]] = i;
		}

		_wantedIdsIndex.resize(_matrixSize);
		for(size_t i=0; i<_matrixSize; i++){
			_wantedIdsIndex[i] = _idToIndex[_wantedIds[i]];
		}
	}

	void writeMatrices(){
		vector<string> matrixFilenames = System::file().listdir(_dirBinaryMatrices);


		for(size_t i=0; i<matrixFilenames.size(); i++){
			string matrixFilename = _dirBinaryMatrices + "/" + matrixFilenames[i];
			if(matrixFilename.find("mat_") == string::npos) continue;

			cout << matrixFilenames[i] << endl;
			writeMatrixASCII(matrixFilenames[i], matrixFilename);
		}

	}

	void writeMatrixASCII(const string& distanceName, const string& binaryMatrixFilename){

		vector<float> rowData(_matrixSize, 0);
		ifstream binaryMatrixFile(binaryMatrixFilename.c_str(), ios::binary);
		string filename = _outputDir + "/" + distanceName + ".csv";
		gzFile out = gzopen((filename + ".gz").c_str(),"wb");

		string str = "";

		for(size_t i=0; i<_matrixSize; i++){
			str += ";" + _wantedIds[i];
		}
		str += '\n';
		gzwrite(out, str.c_str(), str.size());


		for(size_t i=0; i<_matrixSize; i++){

			str = "";
			str += _wantedIds[i] + ";";

			size_t rowIndex = _wantedIdsIndex[i];
			SimkaDistanceMatrixBinary::loadRow(rowIndex, binaryMatrixFile, rowData);

			for(size_t j=0; j<_matrixSize; j++){

				str += Stringify::format("%f", rowData[_wantedIdsIndex[j]]) + ";";

			}

			str.erase(str.size()-1);
			str += '\n';

			gzwrite(out, str.c_str(), str.size());
		}


		gzclose(out);
		binaryMatrixFile.close();
	}


	/*
	void outputMatrixToAscii(){
		string filename = outputDir + "/" + outputFilename + _outputFilenameSuffix + ".csv";
		//IFile* file = gatb::core::system::impl::System::file().newFile(filename, "wb");
		Bag<string>* bag = new BagGzFile<string>(outputDir + "/" + outputFilename + _outputFilenameSuffix + ".csv.gz");
		//Bag<string>* cachedBag = new BagCache<string>(bag, 10000);


		gzFile out = gzopen((filename + ".gz").c_str(),"wb");
		//char buf[BUFSIZ] = { 0 };

		//char buffer[200];
		string str;

		for(size_t i=0; i<matrix.size(); i++){
			str += ";" + bankNames[i];
			//str += ";" + datasetInfos[i]._name;
		}
		str += '\n';
		gzwrite(out, str.c_str(), str.size());
		//cachedBag->insert(str);
		//file->fwrite(str.c_str(), str.size(), 1);

		str = "";

		for(size_t i=0; i<matrix.size(); i++){

			str = "";
			str += bankNames[i] + ";";
			//str += datasetInfos[i]._name + ";";
			for(size_t j=0; j<matrix.size(); j++){

				//snprintf(buffer,200,"%.2f", matrix[i][j]);
				//snprintf(buffer,200,"%f", matrix[i][j]);
				str += Stringify::format("%f", matrix[i][j]) + ";";

				//str += to_string(matrix[i][j]) + ";";
			}

			//matrixNormalizedStr.erase(matrixNormalizedStr.end()-1);
			str.erase(str.size()-1);
			//str.pop_back(); //remove ; at the end of the line
			str += '\n';

			//cout << str << endl;
			//file->fwrite(str.c_str(), str.size(), 1);
			gzwrite(out, str.c_str(), str.size());
			//cachedBag->insert(str);
		}


		//gatb::core::system::IFile* file = gatb::core::system::impl::System::file().newFile(outputDir + "/" + outputFilename + _outputFilenameSuffix + ".csv", "wb");

		//file->flush();
		//delete file;
		gzclose(out);
		//delete cachedBag;

	}*/


};











class SimkaDistanceExporter: public Tool{
public:

	string _execFilename;

	SimkaDistanceExporter(string execFilename): Tool ("Simka-DistanceMatrixWriter"){
		_execFilename = execFilename;

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    parser->push_front (new OptionOneParam (STR_URI_INPUT, "path to a dir of binary matrices", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_INPUT_IDS, "filename of ids in the result matrix (one id per line). Do not used this option to used all ids.", false));
	    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output dir for distance matrices", true));

	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for merged kmer spectrum", true));
	   // parser->push_front (new OptionOneParam (STR_SIMKA2_DATABASE_DIR, "dir path to a simka database", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_NB_PARTITION, "number of partitions", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_PARTITION_ID, "number of the partition", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_MATRIX_DIR, "input filename of k-mer spectrums for which distances has to be computed", true));

	   // IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    //kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", true));

	    //parser->getParser(STR_NB_CORES)->setVisible(false);
		//parser->push_back(kmerParser);

	}



	void execute ()
	{
		IProperties* input = getInput();
		//Parameter params(*this, getInput());
		//Parameter params(input, _execFilename);

		//size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

	    //Integer::apply<Functor,Parameter> (kmerSize, params);

		SimkaDistanceExporterAlgorithm* algo = new SimkaDistanceExporterAlgorithm(input);
		algo->execute();
		delete algo;
	}


};










int main (int argc, char* argv[])
{
    try
    {
    	SimkaDistanceExporter(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



