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

#include "SimkaPotara.hpp"

/*
TODO:
	- Faire la config dans un job a part (job_count.bash) pour avoir la même config pour les job de comptage et la config
	- Verifier les paramètre passer au jobs généré (nbcores, maxmemory...)

*/
SimkaPotara::SimkaPotara(const string& execFilename)  : Tool ("Simka")
{

	_execFilename = execFilename;

	Simka::createOptionsParser(getParser());

	//Kmer parser
    IOptionsParser* coreParser = getParser()->getParser("core");

    //clusterParser->push_back (new OptionNoParam (STR_SIMKA_CLUSTER_MODE, "enable cluster mode. All cluster args below must be set", false));
    coreParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_COUNT, "maximum number of simultaneous counting jobs (a higher value improve execution time but increase temporary disk usage)", false));
    coreParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_MERGE, "maximum number of simultaneous merging jobs (1 job = 1 core)", false));


    IOptionsParser* clusterParser = new OptionsParser ("cluster");
    //clusterParser->push_back (new OptionOneParam (STR_SIMKA_NB_PARTITIONS, "nb partitions", false, "0" ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_COMMAND, "command to submit counting job", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_COMMAND, "command to submit merging job", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_FILENAME, "filename to the couting job template", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_FILENAME, "filename to the merging job template", false ));


	//getParser()->push_back(coreParser);
	getParser()->push_back(clusterParser);
	//getParser()->getParser("core")->getParser(STR_NB_CORES)->setHelp("number of cores per counting job");
    //if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_MAX_MEMORY)))  {  p->setHelp("max memory per counting job (in MBytes) "); }

    //if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_NB_CORES)))  {  p->setHelp("number of cores per job"); }
    //coreParser->push_back(new OptionOneParam(parser->getParser(STR_NB_CORES)->getName(), parser->getParser(STR_NB_CORES)->getHelp(), false, "0"));

	//if (IOptionsParser* input = dskParser->getParser (STR_KMER_ABUNDANCE_MIN_THRESHOLD))  {  input->setVisible (false);  }

	/*
	IOptionsParser* parser = getParser();
	IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();
	parser->push_back (dskParser, 1);
	parser->push_back(dskParser);


	parser->getParser (STR_URI_INPUT)->setHelp("input file of datasets and their id. One dataset per line: dataset_id dataset_filename");
	parser->getParser (STR_KMER_ABUNDANCE_MIN_THRESHOLD)->setVisible (false);
	parser->getParser (STR_HISTOGRAM_MAX)->setVisible (false);
	parser->getParser (STR_URI_SOLID_KMERS)->setVisible (false);
	parser->getParser (STR_URI_OUTPUT_DIR)->setHelp("output directory for temporary files");
	parser->getParser (STR_URI_OUTPUT)->setHelp("output directory for result files");
	parser->getParser (STR_SOLIDITY_KIND)->setHelp("TODO");
	parser->getParser (STR_MINIMIZER_TYPE)->setVisible (false);
	parser->getParser (STR_MINIMIZER_SIZE)->setVisible (false);
	parser->getParser (STR_REPARTITION_TYPE)->setVisible (false);
    if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }

    parser->push_back (new OptionNoParam (STR_SOLIDITY_PER_DATASET.c_str(), "Do not take into consideration multi-counting when determining solidity of kmers", false ));
	*/
	/*
    parser->push_back (new OptionOneParam (STR_URI_INPUT,         "reads file", true ));
    parser->push_back (new OptionOneParam (STR_KMER_SIZE,         "size of a kmer",                           false,  "31"    ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN,"min abundance threshold for solid kmers",  false,  "3"     ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX,"min abundance threshold for solid kmers",  false,  "3"     ));
    parser->push_back (new OptionOneParam (STR_MAX_MEMORY,        "max memory (in MBytes)",                   false, "2000"));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT_DIR,   "output folder for solid kmers",              false));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT,        "output file",                              false));
	*/

    //setParser (parser);
}

struct Parameter
{
    //Parameter (Simka& simka, IProperties* props) : props(props) {}
    Parameter (IProperties* props, const string& execFilename) : _props(props), _execFilename(execFilename) {}
    //Simka& _simka;
    IProperties* _props;
    string _execFilename;
    /*
    string _inputFilename;
    string _outputDir;
    size_t _kmerSize;
    pair<size_t, size_t> _abundanceThreshold;
    bool _soliditySingle;*/
};

template<size_t span> struct Functor  {  void operator ()  (Parameter p)
{
	/*
	cout << "SimkaAlgo.cpp 1" << endl;
	clear();
	delete _banks;
	cout << "SimkaAlgo.cpp 2" << endl;
	SimkaFusion<span>* simkaFusion = new SimkaFusion<span>(_options, _inputFilename, _outputDir, _outputDirTemp, _nbReadsPerDataset, _maxNbReads);
	simkaFusion->execute();
	return;*/

	SimkaPotaraAlgorithm<span> simkaAlgorithm (p._props, p._execFilename);
	simkaAlgorithm.execute();

	/*
#ifdef SIMKA_MIN
	simkaAlgorithm.executeSimkamin();
#else
#endif*/
}};

void SimkaPotara::execute ()
{
	IProperties* input = getInput();
	//Parameter params(*this, getInput());
	Parameter params(input, _execFilename);

	size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    Integer::apply<Functor,Parameter> (kmerSize, params);
}





int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.

		//cout << argv[0] << endl;
        SimkaPotara(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
