/*
 * SimkaMin.cpp
 *
 *  Created on: 24 juin 2017
 *      Author: gbenoit
 */

#include "SimkaMinCount.hpp"
#include "SimkaMinDistance.hpp"
#include "SimkaMinDistanceMatrixExporter.hpp"

void displayHelp(){
	cout << "Usage: ./simkaMin [count|merge|distance|info|export]" << endl;

	cout << endl << "[Distance computation options]" << endl;
	cout << "\tspectrum    : transform datasets in k-mer spectrums" << endl;
	cout << "\tmerge       : merge multiple k-mer spectrum files into a single" << endl;
	cout << "\tdistance    : compute Jaccard and Bray-Curtis distances between k-mer spectrums" << endl;

	cout << endl << "[Distance matrix manipulation options]" << endl;
	cout << "\texport      : export distance matrices stored in binary format" << endl;

	cout << endl << "[Other options]" << endl;
	cout << "\tinfo        : list datasets contained in a k-mer spectrum file" << endl;

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

    		if(programName == "spectrum"){
    			Simka2ComputeKmerSpectrum().run (argc, args);
    		}
    		else if(programName == "merge"){
    			cout << "todo" << endl;
    		}
    		else if(programName == "distance"){
    			SimkaMinDistance().run(argc, args);
    		}
    		else if(programName == "export"){
    			SimkaMinDistanceMatrixExporter().run(argc, args);
    		}
    		else if(programName == "info"){
    			cout << "todo" << endl;
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
