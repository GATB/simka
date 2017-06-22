/*****************************************************************************
 *   Simka: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2015  INRIA
 *   Authors: G.Benoit, C.Lemaitre, P.Peterlongo
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

#ifndef TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_
#define TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_

#include <gatb/gatb_core.hpp>
#include "../export/SimkaDistanceMatrixBinary.hpp"

const string STR_SIMKA_DISTANCE_BRAYCURTIS = "-bray-curtis";
const string STR_SIMKA_DISTANCE_CHORD = "-chord";
const string STR_SIMKA_DISTANCE_HELLINGER = "-hellinger";
const string STR_SIMKA_DISTANCE_CANBERRA = "-canberra";
const string STR_SIMKA_DISTANCE_KULCZYNSKI = "-kulczynski";


#define MULTISCALE_BOOSTRAP_NB_BOOSTRAPS 1

typedef vector<u_int16_t> SpeciesAbundanceVectorType;

enum SIMKA_MATRIX_TYPE{
	SYMETRICAL,
	ASYMETRICAL,
};

class DistanceMatrixData
{

public:

	//size_t _nbBanks;
	//size_t _nbNewBanks;

	vector<vector<u_int64_t> > _matrix_rectangular;
	vector<vector<u_int64_t> > _matrix_squaredHalf;
	vector<u_int64_t> _marginalValues;

	DistanceMatrixData(){};

	void resize(u_int64_t nbBanks, u_int64_t nbNewBanks){

		//cout << "Create distance matrix data" << endl;
		//cout << "\t" << nbBanks << " " << nbNewBanks << endl << endl;
		//_nbBanks = nbBanks;
		//_nbNewBanks = nbNewBanks;

		u_int64_t nbOldBanks = nbBanks - nbNewBanks;

		if(nbOldBanks > 0){
			_matrix_rectangular.resize(nbNewBanks);
		}
		else{
			_matrix_rectangular.resize(0);
		}

		_matrix_squaredHalf.resize(nbNewBanks-1);

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			_matrix_rectangular[i].resize(nbOldBanks);
		}

		u_int64_t size = nbNewBanks-1;
		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			_matrix_squaredHalf[i].resize(size);
			size -= 1;
		}

		_marginalValues.resize(nbNewBanks, 0);
	}

	DistanceMatrixData& operator+=  (const DistanceMatrixData& other){

		for(size_t i=0; i < _marginalValues.size(); i++){
			_marginalValues[i] += other._marginalValues[i];
		}

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				_matrix_rectangular[i][j] += other._matrix_rectangular[i][j];
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				_matrix_squaredHalf[i][j] += other._matrix_squaredHalf[i][j];
			}
		}

		return *this;
	}

	void load(Iterator<long double>* gzIt){

		for(size_t i=0; i < _marginalValues.size(); i++){
			_marginalValues[i] = gzIt->item();
			gzIt->next();
		}

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				_matrix_rectangular[i][j] = gzIt->item();
				gzIt->next();
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				_matrix_squaredHalf[i][j] = gzIt->item();
				gzIt->next();
			}
		}
	}

	void save (Bag<long double>* bag){

		for(size_t i=0; i < _marginalValues.size(); i++){
			bag->insert((long double)_marginalValues[i]);
		}

		for(size_t i=0; i<_matrix_rectangular.size(); i++){
			for(size_t j=0; j<_matrix_rectangular[i].size(); j++){
				bag->insert((long double)_matrix_rectangular[i][j]);
			}
		}

		for(size_t i=0; i<_matrix_squaredHalf.size(); i++){
			for(size_t j=0; j<_matrix_squaredHalf[i].size(); j++){
				bag->insert((long double)_matrix_squaredHalf[i][j]);
			}
		}

	}

};
/*
class SimkaDistanceParam{

public:

	SimkaDistanceParam(){}

	SimkaDistanceParam(IProperties* params){
	    //_computeBrayCurtis = true;
		//_computeChord = true;
		//_computeHellinger = true;
		//_computeCanberra = true;
		//_computeKulczynski = true;
		//_computeBrayCurtis = params->get(STR_SIMKA_DISTANCE_BRAYCURTIS);
		//_computeChord = params->get(STR_SIMKA_DISTANCE_CHORD);
		//_computeHellinger = params->get(STR_SIMKA_DISTANCE_HELLINGER);
		//_computeCanberra = params->get(STR_SIMKA_DISTANCE_CANBERRA);
		//_computeKulczynski = params->get(STR_SIMKA_DISTANCE_KULCZYNSKI);
	}

	//bool _computeBrayCurtis;
	//bool _computeChord;
	//bool _computeHellinger;
	//bool _computeCanberra;
	//bool _computeKulczynski;
};*/


class SimkaStatistics{

public:

	SimkaStatistics(size_t nbBanks, size_t nbNewBanks, bool computeSimpleDistances, bool computeComplexDistances);
	SimkaStatistics& operator+=  (const SimkaStatistics& other);
	void print();
	void load(const string& filename);
	void save(const string& filename);
	void outputMatrix(const string& outputDirTemp, const vector<string>& _bankNames, u_int64_t nbProcessedDatasets, u_int64_t nbNewDatasets);

	void loadInfoSimka1(const string& tmpDir, const vector<string>& datasetIds){
		_totalReads = 0;

		for(size_t i=0; i<_nbBanks; i++){

			string name = datasetIds[i];
			string countFilename = tmpDir + "/count_synchro/" +  name + ".ok";

			string line;
			ifstream file(countFilename.c_str());
			vector<string> lines;
			while(getline(file, line)){
				if(line == "") continue;
				lines.push_back(line);
			}
			file.close();

			u_int64_t nbReads = strtoull(lines[0].c_str(), NULL, 10);


			_datasetNbReads[i] = nbReads;
			_nbSolidDistinctKmersPerBank[i] = strtoull(lines[1].c_str(), NULL, 10);
			_nbSolidKmersPerBank[i] = strtoull(lines[2].c_str(), NULL, 10);


			if(_computeSimpleDistances){
				_chord_sqrt_N2[i] = sqrt(strtoull(lines[3].c_str(), NULL, 10));
			}

			_totalReads += nbReads;
			/*
			for (size_t j=0; j<_nbCores; j++){
				DistanceCommand<span>* cmd = dynamic_cast<DistanceCommand<span>*>(_cmds[j]);
				cmd->_stats->_datasetNbReads[i] = nbReads;
				cmd->_stats->_nbSolidDistinctKmersPerBank[i] = strtoull(lines[1].c_str(), NULL, 10);
				cmd->_stats->_nbSolidKmersPerBank[i] = strtoull(lines[2].c_str(), NULL, 10);
				cmd->_stats->_chord_sqrt_N2[i] = sqrt(strtoull(lines[3].c_str(), NULL, 10));
			}*/
		}
	}

    size_t _nbBanks;
    size_t _nbNewBanks;
    size_t _symetricDistanceMatrixSize;
    bool _computeSimpleDistances;
    bool _computeComplexDistances;

    double _totalReads;

	vector<u_int64_t> _nbSolidDistinctKmersPerBank;
	vector<u_int64_t> _nbSolidKmersPerBank;

	vector<u_int64_t> _nbDistinctKmersSharedByBanksThreshold;
	vector<u_int64_t> _nbKmersSharedByBanksThreshold;

	//vector<vector<u_int64_t> > _matrixNbDistinctSharedKmers;
	vector<vector<u_int64_t> > _matrixNbSharedKmers;

	DistanceMatrixData _brayCurtisNumerator;
	DistanceMatrixData _matrixNbDistinctSharedKmers;


	vector<vector<DistanceMatrixData> > _multiscaleBoostrap_distanceData_braycurtis;
	//vector<vector<u_int64_t> > _brayCurtisNumerator;
	//vector<vector<double> > _kullbackLeibler;


    //Abundance Chord
	vector<vector<long double> > _chord_NiNj;
	vector<long double> _chord_sqrt_N2;

    //Abundance Hellinger
	vector<vector<u_int64_t> > _hellinger_SqrtNiNj;
	vector<vector<u_int64_t> > _whittaker_minNiNj;
	vector<vector<long double> > _kullbackLeibler;
	vector<vector<u_int64_t> > _abundance_jaccard_intersection;

    //Abundance Canberra
	vector<vector<u_int64_t> > _canberra;

	//Abundance Kulczynski
	vector<vector<u_int64_t> > _kulczynski_minNiNj;

	//string _outputDir;

	u_int64_t _nbKmers;
	vector<u_int64_t> _nbKmersPerBank;
	u_int64_t _nbErroneousKmers;

	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSolidKmers;

	u_int64_t _nbSharedKmers;
	u_int64_t _nbDistinctSharedKmers;
	//SimkaDistanceParam _distanceParams;

	vector<u_int64_t> _datasetNbReads;
	//u_int64_t _nbKmersInCoupleBankSupRatio;

	//unordered_map<string, histo_t> _histos;


private:

	string _outputFilenameSuffix;
};


class SimkaDistance {

public:

	vector<vector<float> > _matrix; //to remove
	vector<vector<float> > _matrix_rectangular;
	vector<vector<float> > _matrix_squaredHalf;


	SimkaDistance(SimkaStatistics& stats);

	void clearMatrix(){
		/*
		for(size_t i=0; i<_matrix.size(); i++){
			for(size_t j=0; j<_matrix.size(); j++){
				_matrix[i][j] = 0;
			}
		}*/
	}

    void _matrixJaccardAbundance(){
    	clearMatrix();

    	//for(size_t i=0; i<_nbBanks; i++){
    	//	for(size_t j=i+1; j<_nbBanks; j++){
		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_jaccard(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixBrayCurtis(size_t ii, size_t jj){

    	clearMatrix();

		size_t nbOldBanks = _nbBanks - _nbNewBanks;

		//cout << "lala" << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular[0].size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf[0].size() << endl;

		for(size_t i=0; i<_stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_rectangular[i].size(); j++){

				double dist = distance_abundance_brayCurtis(i+nbOldBanks, j, i, j, _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_rectangular, _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._marginalValues);
				_matrix_rectangular[i][j] = dist;
			}
		}

		for(size_t i=0; i<_stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_squaredHalf.size(); i++){

			u_int64_t jOffset = _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_squaredHalf.size() - _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_squaredHalf[i].size();

			for(size_t j=0; j<_stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_squaredHalf[i].size(); j++){

				double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._matrix_squaredHalf, _stats._multiscaleBoostrap_distanceData_braycurtis[ii][jj]._marginalValues);
				_matrix_squaredHalf[i][j] = dist;
			}
		}

		/*
		cout << endl;
		cout << endl;
		for(size_t i=0; i<_stats._nbSolidKmersPerBank.size(); i++){
			cout << _stats._nbSolidKmersPerBank[i] << " ";
		}
		cout << endl;
		cout << endl;
		for(size_t i=0; i<_stats._brayCurtisNumerator._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._brayCurtisNumerator._matrix_rectangular[i].size(); j++){
				cout << _stats._brayCurtisNumerator._matrix_rectangular[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;


		cout << endl;
		for(size_t i=0; i<_stats._brayCurtisNumerator._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._brayCurtisNumerator._matrix_rectangular[i].size(); j++){
				cout << _matrix_rectangular[i][j] << " ";
			}
			cout << endl;
		}*/

		//cout << _matrix.size() << " " << _matrix[0].size() << "   " << _nbNewBanks << endl;
		/*
		for(size_t i=0; i<nbOldBanks; i++){
			for(size_t j=0; j<_nbNewBanks; j++){

				double dist = distance_abundance_brayCurtis(i, j, i, j+nbOldBanks, _stats._brayCurtisNumerator._matrix_rectangular);
				_matrix[i][j] = dist;

			}
		}

		cout << "Final matrix size: " << _matrix.size() << " " << _matrix[0].size() << endl;


		u_int64_t jOffset = 1;

		for(size_t i=nbOldBanks; i<_nbBanks; i++){
			cout << i << endl;
			for(size_t j=(i-nbOldBanks)+1; j<_nbNewBanks; j++){
				cout << "\t" << j << endl;

				double dist = distance_abundance_brayCurtis(i-nbOldBanks-jOffset, j, i, j+nbOldBanks, _stats._brayCurtisNumerator._matrix_squaredHalf);
				_matrix[i][j] = dist;


				size_t iTemp = i-nbOldBanks;
				size_t i2 = j;
				size_t j2 = iTemp;
				i2 += nbOldBanks;
				_matrix[i2][j2] = dist;

			}

			jOffset += 1;
		}
		*/

		/*
		u_int64_t jOffset = 1;

		for(size_t i=0; i<_nbNewBanks; i++){
			cout << i << endl;
			for(size_t j=i+1; j<_nbNewBanks; j++){
				cout << "\t" << j << endl;

				double dist = distance_abundance_brayCurtis(i+nbOldBanks, j-jOffset, i+nbOldBanks, j, _stats._brayCurtisNumerator._matrix_squaredHalf);
				_matrix[i+nbOldBanks][j] = dist;


				size_t iTemp = i-nbOldBanks;
				size_t i2 = j;
				size_t j2 = iTemp;
				i2 += nbOldBanks;
				_matrix[i2][j2] = dist;

			}

			jOffset += 1;
		}*/

    }

    void _matrixChord(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_chord(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixHellinger(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_hellinger(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixWhittaker(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_whittaker(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixKullbackLeibler(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_kullbackLeibler(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixCanberra(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			get_abc(i, j, j, a, b ,c);
    			double dist = distance_abundance_canberra(i, j, a, b, c);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}
		*/

    }

    void _matrixKulczynski(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_kulczynski(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixSymJaccardAbundance(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_jaccard_simka(i, j, SYMETRICAL);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixAsymJaccardAbundance(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			_matrix[i][j] = distance_abundance_jaccard_simka(i, j, ASYMETRICAL);
    			_matrix[j][i] = distance_abundance_jaccard_simka(j, i, ASYMETRICAL);
    		}
    	}


    }

    void _matrixOchiai(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_ochiai(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrixSorensen(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double dist = distance_abundance_sorensen(i, j);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrix_presenceAbsence_sorensenBrayCurtis(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			get_abc(i, j, j, a, b ,c);

    			double dist = distance_presenceAbsence_sorensenBrayCurtis(a, b, c);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}
		*/

    }

    void _matrix_presenceAbsence_Whittaker(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			get_abc(i, j, j, a, b ,c);

    			double dist = distance_presenceAbsence_whittaker(a, b, c);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}
		*/

    }

    void _matrix_presenceAbsence_kulczynski(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			get_abc(i, j, j, a, b ,c);

    			double dist = distance_presenceAbsence_kulczynski(a, b, c);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}
		*/

    }

    void _matrix_presenceAbsence_ochiai(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			get_abc(i, j, j, a, b ,c);

    			double dist = distance_presenceAbsence_ochiai(a, b, c);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}
		*/

    }

    void _matrix_presenceAbsence_chordHellinger(){
    	/*
    	clearMatrix();
    	u_int64_t a, b, c;

    	if(_nbBanks == _nbNewBanks){

			for(size_t i=0; i<_nbBanks; i++){
				for(size_t j=i+1; j<_nbBanks; j++){

        			get_abc(i, j, j, a, b ,c);

        			double dist = distance_presenceAbsence_chordHellinger(a, b, c);
        			//cout << a << "  " << b << " " << c << "  " << dist << endl;
        			_matrix[i][j] = dist;
        			_matrix[j][i] = dist;
        		}
        	}

    	}
    	else{

    		size_t nbOldBanks = _nbBanks - _nbNewBanks;

    		//cout << _matrix.size() << " " << _matrix[0].size() << "   " << _nbNewBanks << endl;
			for(size_t i=0; i<nbOldBanks; i++){
				for(size_t j=0; j<_nbNewBanks; j++){

    				//cout << i << " " << j << endl;
        			get_abc(i, j, j+nbOldBanks, a, b ,c);

        			double dist = distance_presenceAbsence_chordHellinger(a, b, c);
        			_matrix[i][j] = dist;

        		}
        	}

			for(size_t i=nbOldBanks; i<_nbBanks; i++){
				for(size_t j=(i-nbOldBanks)+1; j<_nbNewBanks; j++){

        			get_abc(i, j, j+nbOldBanks, a, b ,c);

        			double dist = distance_presenceAbsence_chordHellinger(a, b, c);
        			_matrix[i][j] = dist;


        			size_t iTemp = i-nbOldBanks;
        			size_t i2 = j;
        			size_t j2 = iTemp;
        			i2 += nbOldBanks;
        			_matrix[i2][j2] = dist;

        		}
        	}
    	}
		*/
    }

    void _matrix_presenceAbsence_jaccardCanberra(){
    	u_int64_t a, b, c;

    	clearMatrix();

		size_t nbOldBanks = _nbBanks - _nbNewBanks;

		//cout << "lala" << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_rectangular[0].size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf.size() << endl;
		//cout << _stats._brayCurtisNumerator._matrix_squaredHalf[0].size() << endl;

		for(size_t i=0; i<_stats._matrixNbDistinctSharedKmers._matrix_rectangular.size(); i++){
			for(size_t j=0; j<_stats._matrixNbDistinctSharedKmers._matrix_rectangular[i].size(); j++){


				get_abc(i+nbOldBanks, j, i, j, _stats._matrixNbDistinctSharedKmers._matrix_rectangular, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);

				//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j, i, j, _stats._brayCurtisNumerator._matrix_rectangular);
				_matrix_rectangular[i][j] = dist;
			}
		}

		for(size_t i=0; i<_stats._matrixNbDistinctSharedKmers._matrix_squaredHalf.size(); i++){

			u_int64_t jOffset = _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf.size() - _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf[i].size();

			for(size_t j=0; j<_stats._matrixNbDistinctSharedKmers._matrix_squaredHalf[i].size(); j++){

				get_abc(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._matrixNbDistinctSharedKmers._matrix_squaredHalf, a, b ,c);
				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);
				//double dist = distance_abundance_brayCurtis(i+nbOldBanks, j+nbOldBanks+1+jOffset, i, j, _stats._brayCurtisNumerator._matrix_squaredHalf);
				_matrix_squaredHalf[i][j] = dist;
			}
		}

		/*
		//cout << _matrix.size() << " " << _matrix[0].size() << "   " << _nbNewBanks << endl;
		for(size_t i=0; i<nbOldBanks; i++){
			for(size_t j=0; j<_nbNewBanks; j++){

				//cout << i << " " << j << endl;
				get_abc(i, j, j+nbOldBanks, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);
				_matrix[i][j] = dist;

			}
		}

		for(size_t i=nbOldBanks; i<_nbBanks; i++){
			for(size_t j=(i-nbOldBanks)+1; j<_nbNewBanks; j++){

				get_abc(i, j, j+nbOldBanks, a, b ,c);

				double dist = distance_presenceAbsence_jaccardCanberra(a, b, c);
				_matrix[i][j] = dist;


				size_t iTemp = i-nbOldBanks;
				size_t i2 = j;
				size_t j2 = iTemp;
				i2 += nbOldBanks;
				_matrix[i2][j2] = dist;

			}
		}*/



    }

    void _matrix_presenceAbsence_jaccard_simka(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			double dist = distance_presenceAbsence_jaccard_simka(i, j, SYMETRICAL);
    			_matrix[i][j] = dist;
    			_matrix[j][i] = dist;
    		}
    	}


    }

    void _matrix_presenceAbsence_jaccard_simka_asym(){
    	clearMatrix();

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){

    			_matrix[i][j] = distance_presenceAbsence_jaccard_simka(i, j, ASYMETRICAL);
    			_matrix[j][i] = distance_presenceAbsence_jaccard_simka(j, i, ASYMETRICAL);
    		}
    	}


    }


    void computeJaccardDistanceFromBrayCurtis(const vector<vector<float> >& brayDistanceMatrix){
    	//vector<vector<float> > jaccardDistanceMatrix = createSquaredMatrix(_nbBanks, _nbNewBanks);

		for(size_t j=0; j<_nbNewBanks; j++){
			for(size_t i=0; i<j; i++){
    			double B = brayDistanceMatrix[i][j];
    			double J = (2*B) / (1+B);
    			_matrix[i][j] = J;
    		}
    	}

    }

    /*
    static void SimkaDistance::createSquaredMatrix(size_t n, size_t m, vector<vector<float> >& matrix){

        matrix.resize(n);
        for(size_t i=0; i<n; i++)
        	matrix[i].resize(m, 0);

        //

    }*/

private:


	void get_abc(size_t bank1, size_t bank2, size_t i2, size_t j2, vector<vector<u_int64_t> >& crossedData, u_int64_t& a, u_int64_t& b, u_int64_t& c);


    double distance_abundance_brayCurtis(size_t i, size_t j, size_t i2, size_t j2, vector<vector<u_int64_t> >& crossedData, vector<u_int64_t> marginalData){

    	//double intersection = _stats._abundance_jaccard_intersection[i][j];
    	double union_ = marginalData[i] + marginalData[j];

    	if(union_ == 0) return 1;

    	double intersection = 2 * crossedData[i2][j2];

    	double jaccard = 1 - intersection / union_;

    	//cout << jaccard << "   " << _stats._nbSolidKmersPerBank[i] << " " <<  _stats._nbSolidKmersPerBank[j] << "  " << _stats._brayCurtisNumerator[i][j] << endl;
    	return jaccard;
    }

    /*
    double distance_abundance_brayCurtis_squaredHalf(size_t i, size_t j, size_t j2){

    	//double intersection = _stats._abundance_jaccard_intersection[i][j];
    	double union_ = _stats._nbSolidKmersPerBank[i] + _stats._nbSolidKmersPerBank[j2];

    	if(union_ == 0) return 1;

    	double intersection = 2 * _stats._brayCurtisNumerator[i][j];

    	double jaccard = 1 - intersection / union_;

    	//cout << jaccard << "   " << _stats._nbSolidKmersPerBank[i] << " " <<  _stats._nbSolidKmersPerBank[j] << "  " << _stats._brayCurtisNumerator[i][j] << endl;
    	return jaccard;
    }*/


    double distance_abundance_chord(size_t i, size_t j);
    double distance_abundance_hellinger(size_t i, size_t j);
    //double distance_abundance_jaccard_intersection(size_t i, size_t j);
    double distance_abundance_whittaker(size_t i, size_t j);
    double distance_abundance_kullbackLeibler(size_t i, size_t j);
    double distance_abundance_canberra(size_t i, size_t j, u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_abundance_kulczynski(size_t i, size_t j);
    double distance_abundance_ochiai(size_t i, size_t j);
    double distance_abundance_sorensen(size_t i, size_t j);
    double distance_abundance_jaccard(size_t i, size_t j);
    double distance_abundance_jaccard_simka(size_t i, size_t j, SIMKA_MATRIX_TYPE type);


    double distance_presenceAbsence_chordHellinger(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_hellinger(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_whittaker(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_canberra(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_kulczynski(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_ochiai(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_sorensenBrayCurtis(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_jaccardCanberra(u_int64_t& ua, u_int64_t& ub, u_int64_t& uc);
    double distance_presenceAbsence_jaccard_simka(size_t i, size_t j, SIMKA_MATRIX_TYPE type);

	SimkaStatistics& _stats;
	//SimkaDistanceParam _distanceParams;
	size_t _nbBanks;
	size_t _nbNewBanks;

};




#endif /* TOOLS_SIMKA_SRC_SIMKADISTANCE_HPP_ */
