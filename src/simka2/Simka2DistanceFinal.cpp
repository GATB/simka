/*
 * Simka2ComputeKmerSpectrum.hpp
 *
 *  Created on: 4 nov. 2016
 *      Author: gbenoit
 */

#include <gatb/gatb_core.hpp>
//#include "../core/SimkaUtils.hpp"
#include "Simka2Utils.hpp"
#include "../core/SimkaUtils.hpp"
#include "Simka2Database.hpp"
#include "../export/SimkaDistanceMatrixBinary.hpp"
//#include "../minikc/MiniKC.hpp"
//#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
//#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
//#include "SimkaAlgorithm.hpp"
//#include "SimkaAlgorithm.hpp"




/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class Simka2DistanceFinalAlgorithm : public Algorithm
{

public:

	string _databaseDir;
	//string _inputFilenameUnknown;
	string _outputDir;
	size_t _nbCores;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	size_t _kmerSize;

	Simka2Database _database;
	size_t _nbPartitions;

	SimkaStatistics* _stats;

	size_t _nbBanks;
	size_t _nbNewBanks;
	size_t _nbOldBanks;

	string _dirMatrixParts;
	string _dirMatrixMatrixBinary;
	string _dirMatrixMatrixBinaryTmp;

	vector<vector<float> > _matrix;
	vector<string> _oldIds;
	vector<string> _newIds;

	Simka2DistanceFinalAlgorithm(IProperties* options):
		Algorithm("simka", -1, options)
	{

	}

	void execute(){
		parseArgs();
		createDatabases();
		//initStatistics();
		distance();

		cout << _oldIds.size() << " " << _newIds.size() << endl;
		_oldIds.insert( _oldIds.end(), _newIds.begin(), _newIds.end() ); //concat two vectors
		cout << _oldIds.size() << " " << _newIds.size() << endl;
		SimkaDistanceMatrixBinary::saveMatrixIds(_dirMatrixMatrixBinary, _oldIds, _oldIds.size());
	}

	void parseArgs(){

    	_nbCores =  getInput()->getInt(STR_NB_CORES);
    	_kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	//_partitionId = getInput()->getInt(STR_SIMKA2_PARTITION_ID);
    	_databaseDir =  getInput()->getStr(STR_SIMKA2_DATABASE_DIR);
    	_nbPartitions =  getInput()->getInt(STR_SIMKA2_NB_PARTITION);


    	_dirMatrixParts = _databaseDir + "/distance/temp_parts";
    	_dirMatrixMatrixBinary = _databaseDir + "/distance/matrix_binary";
    	_dirMatrixMatrixBinaryTmp = _databaseDir + "/distance/matrix_binary_temp";
	}

	void createDatabases(){

		_database = Simka2Database(_databaseDir);

		_nbBanks = _database._entries.size();
		_nbNewBanks = _nbBanks - _database._nbProcessedDataset;

	}

	void distance(){
		cout << endl << "Computing stats..." << endl;
		//cout << this->_nbBanks << endl;

		//u_int64_t nbKmers = 0;

		//SimkaDistanceParam distanceParams(this->_options);
		vector<string> kmerSpectrumDirs;
		SimkaStatistics mainStats(this->_nbBanks, this->_nbNewBanks, this->_computeSimpleDistances, this->_computeComplexDistances);
		simka2_loadStatInfos(_databaseDir, _database._uniqKmerSpectrumDirs, _database._entries, kmerSpectrumDirs, &mainStats, _database._entriesInfos);

		for(size_t i=0; i<_nbPartitions; i++){

			//cout << mainStats._nbDistinctKmers << endl;
			string filename = _dirMatrixParts + "/" + Stringify::format("%i", i) + ".gz";
			//string filename =  + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".gz";
			//Storage* storage = StorageFactory(STORAGE_HDF5).load (this->_outputDirTemp + "/stats/part_" + SimkaAlgorithm<>::toString(i) + ".stats");
			//LOCAL (storage);

			SimkaStatistics stats(this->_nbBanks, this->_nbNewBanks, this->_computeSimpleDistances, this->_computeComplexDistances);
			kmerSpectrumDirs.clear();
			simka2_loadStatInfos(_databaseDir, _database._uniqKmerSpectrumDirs, _database._entries, kmerSpectrumDirs, &stats, _database._entriesInfos);
			stats.load(filename);

			//cout << stats._nbDistinctKmers << "   " << stats._nbKmers << endl;
			mainStats += stats;

			//nbKmers += stats._nbKmers;
		}


		//cout << mainStats._nbDistinctKmers << endl;
		//cout << "Nb kmers: " << nbKmers << endl;

		//getCountInfo(mainStats);
		//for(size_t i=0; i<this->_nbBanks; i++){
		//	cout << mainStats._nbSolidDistinctKmersPerBank[i] << endl;
		//}
		mainStats.outputMatrix(_dirMatrixMatrixBinaryTmp, _database._entries);

		mergeDistanceMatrix();

		mainStats.print();
	}

	void mergeDistanceMatrix(){

		SimkaDistanceMatrixBinary::loadMatrixIds(_dirMatrixMatrixBinary, _oldIds);
		SimkaDistanceMatrixBinary::loadMatrixIds(_dirMatrixMatrixBinaryTmp, _newIds);

		_nbBanks = _oldIds.size() + _newIds.size();
		_nbOldBanks = _oldIds.size();
		_nbNewBanks = _newIds.size();

		cout << _nbBanks << " " << _nbOldBanks << " " << _nbNewBanks << endl;

		_matrix.resize(_nbBanks);
	    for(size_t i=0; i<_matrix.size(); i++)
	    	_matrix[i].resize(_nbNewBanks, 0);






		vector<string> matrixFilenames = System::file().listdir(_dirMatrixMatrixBinaryTmp);

		//No existing run of simka, we just copy binary matrix
		if(_oldIds.size() == 0){

			for(size_t i=0; i<matrixFilenames.size(); i++){
				string matrixFilenameTemp = _dirMatrixMatrixBinaryTmp + "/" + matrixFilenames[i];
				if(matrixFilenameTemp.find("mat_") == string::npos) continue;

				//cout << matrixFilenames[i] << endl;
				string matrixFilenameExisting = _dirMatrixMatrixBinary + "/" + matrixFilenames[i];


				System::file().rename(matrixFilenameTemp, matrixFilenameExisting);
			}

			System::file().rename(_dirMatrixMatrixBinaryTmp + "/matrix_infos.bin", _dirMatrixMatrixBinary + "/matrix_infos.bin");
		}
		else{


			//vector<string> newMatrixIds;
			//SimkaDistanceMatrixBinary::loadMatrixIds(_dirMatrixMatrixBinaryTmp, newMatrixIds);



			//cout << newMatrixIds.size() << "   " << oldMatrixIds.size() << endl;
			//cout << nbBanks << "  " << nbNewBanks << endl;
			//cout << "ha!" << endl;
			//exit(1);
			for(size_t i=0; i<matrixFilenames.size(); i++){
				string matrixFilenameTemp = _dirMatrixMatrixBinaryTmp + "/" + matrixFilenames[i];
				if(matrixFilenameTemp.find("mat_") == string::npos) continue;

				cout << matrixFilenames[i] << endl;
				string matrixFilenameExisting = _dirMatrixMatrixBinary + "/" + matrixFilenames[i];

				mergeMatrix(matrixFilenames[i], matrixFilenameExisting, matrixFilenameTemp);
			}

		}

	}



	void mergeMatrix(const string& distanceName, const string& oldMatrixFilename, const string& newMatrixFilename){

		//exit(1);
		cout << "merging matrix: " << distanceName << endl;


		ifstream oldMatrixFile(oldMatrixFilename.c_str(), ios::binary);
		ifstream newMatrixFile(newMatrixFilename.c_str(), ios::binary);

		ofstream mergedMatrixFile((oldMatrixFilename + ".temp").c_str(), ios::binary);

		SimkaDistanceMatrixBinary::loadMatrix(newMatrixFilename, _matrix);

		for(size_t i=0; i<_matrix.size(); i++){
			for(size_t j=0; j<_matrix[i].size(); j++){
				cout << _matrix[i][j] << endl;
			}
		}

		u_int64_t oldRowSize = sizeof(float)*_nbOldBanks;
		u_int64_t newRowSize = sizeof(float)*_nbNewBanks;
		u_int64_t rowSize = sizeof(float)*_nbBanks;

		vector<float> oldRowData(_nbOldBanks);
		vector<float> newRowData(_nbNewBanks);
		vector<float> rowData(_nbBanks);

		for(size_t i=0; i<_nbOldBanks; i++){

			oldMatrixFile.read((char*)oldRowData.data(), oldRowSize);
			mergedMatrixFile.write((char*)oldRowData.data(), oldRowSize);

			for(size_t j=0; j<_nbNewBanks; j++){
				newRowData[j] = _matrix[i][j];
			}

			mergedMatrixFile.write((char*)newRowData.data(), newRowSize);
		}

		for(size_t i=0; i<_nbNewBanks; i++){

			//oldMatrixFile.read((char*)oldRowData.data(), oldRowSize);
			//mergedMatrixFile.write((char*)oldRowData.data(), oldRowSize);

			for(size_t j=0; j<_nbBanks; j++){
				rowData[j] = _matrix[j][i];
			}

			mergedMatrixFile.write((char*)rowData.data(), rowSize);
		}

		oldMatrixFile.close();
		newMatrixFile.close();
		mergedMatrixFile.close();

		System::file().remove(oldMatrixFilename);
		System::file().rename(oldMatrixFilename + ".temp", oldMatrixFilename);

	}

	void writeFinishSignal(){

		//string finishFilename = _outputDir + "/" + _datasetID + "-success";
		//IFile* file = System::file().newFile(finishFilename, "w");
		//delete file;
	}



};











class Simka2DistanceFinal: public Tool{
public:

	string _execFilename;

	Simka2DistanceFinal(string execFilename): Tool ("Simka2-Distance"){
		_execFilename = execFilename;

	    IOptionsParser* parser = getParser();//new OptionsParser ("Simka2 - Compute Kmer Spectrum");

	    //parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for merged kmer spectrum", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_DATABASE_DIR, "dir path to a simka database", true));
	    parser->push_front (new OptionOneParam (STR_SIMKA2_NB_PARTITION, "number of partitions", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_PARTITION_ID, "number of the partition", true));
	    //parser->push_front (new OptionOneParam (STR_SIMKA2_DISTANCE_MATRIX_DIR, "input filename of k-mer spectrums for which distances has to be computed", true));

	    IOptionsParser* kmerParser = new OptionsParser ("kmer");
	    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", true));

	    //parser->getParser(STR_NB_CORES)->setVisible(false);
		parser->push_back(kmerParser);

	}


	struct Parameter
	{
	    //Parameter (Simka& simka, IProperties* props) : props(props) {}
	    Parameter (IProperties* props, const string& execFilename) : _props(props), _execFilename(execFilename) {}
	    //Simka& _simka;
	    IProperties* _props;
	    string _execFilename;
	};

	template<size_t span> struct Functor {

    	void operator ()  (Parameter p){

    		Simka2DistanceFinalAlgorithm<span>* algo = new Simka2DistanceFinalAlgorithm<span>(p._props);
    		algo->execute();
    		delete algo;
    	}


	};

	void execute ()
	{
		IProperties* input = getInput();
		//Parameter params(*this, getInput());
		Parameter params(input, _execFilename);

		size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

	    Integer::apply<Functor,Parameter> (kmerSize, params);
	}

};









int main (int argc, char* argv[])
{
    try
    {
    	Simka2DistanceFinal(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



