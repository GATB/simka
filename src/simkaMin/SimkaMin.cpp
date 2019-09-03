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

#include "SimkaMinCount.hpp"
#include "SimkaMinDistance.hpp"
#include "SimkaMinDistanceMatrixExporter.hpp"
#include "SimkaMinDistanceMatrixMerger.hpp"
#include "SimkaMinInfos.hpp"
#include "SimkaMinAppend.hpp"

void displayHelp(){
	cout << "Usage: ./simkaMin [option]" << endl;

	cout << endl << "[Distance computation options]" << endl;
	cout << "\tsketch      : transform datasets in small sketches of k-mers and their abundance" << endl;
	cout << "\tdistance    : compute Jaccard and Bray-Curtis distances between sketches" << endl;

	cout << endl << "[Distance matrix manipulation options]" << endl;
	cout << "\texport      : export distance matrices stored in binary format" << endl;
	//cout << "\tmatrix-update       : update existing distance matrices" << endl;

	cout << endl << "[Sketch options]" << endl;
	cout << "\tappend      : merge multiple sketch files into a single one" << endl;
	cout << "\tinfo        : list datasets contained in a sketch file" << endl;

	cout << endl;
}

int main (int argc, char* argv[])
{
    try
    {
    	if(argc < 2){
    		displayHelp();
    	}
    	else{

    		//std::vector<char*>  args;
    		vector<char*> argsTemp( argv, argv + argc );
    		argsTemp.erase(argsTemp.begin()+1);
    		//std::transform(argsTemp.begin(), argsTemp.end(), std::back_inserter(vc), convert);

    		char** args = &argsTemp[0];
    		//char* args[];

    		//for(string& arg: argsTemp){

    		//}
    		//rArray = new char*[argc+1];
    		//for(int i=0; i <= argc; i++) {
    		//    rArray[i] = argv[i];
    		//}
    		// use rArray
    		//delete [] rArray;


    		//char* args = new char*[argc-1];
    		//vector<string> test;

    		//for(size_t i=0; i<argc; i++){
    		//	if (i==1) continue;
    		//	args[i] = argv[i];
    		//}

    		argc -= 1;
    		//vector<char*> args(argv);

    		string programName = string(argv[1]);

    		if(programName == "sketch"){
    			Simka2ComputeKmerSpectrum().run (argc, args);
    		}
    		else if(programName == "append"){
    			SimkaMinAppend().run(argc, args);
    		}
    		else if(programName == "distance"){
    			SimkaMinDistance().run(argc, args);
    		}
    		else if(programName == "export"){
    			SimkaMinDistanceMatrixExporter().run(argc, args);
    		}
    		else if(programName == "matrix-update"){ //Hidden feature
    			SimkaMinDistanceMatrixMerger().run(argc, args);
    		}
    		else if(programName == "info"){
    			SimkaMinInfos().run(argc, args);
    		}
    		else{
    			displayHelp();
    		}
    	}
    	//cout << argc << endl;
    	//cout << argv[0] << endl;
    	//cout << argv[1] << endl;
    	//cout << argv[2] << endl;
    	//
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
