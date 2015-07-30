
#include "SimkaPotala.hpp"


SimkaPotala::SimkaPotala()  : Tool ("SimkaPotala")
{

	Simka::createOptionsParser(getParser());

	//Kmer parser
    IOptionsParser* jobParser = new OptionsParser ("job");

    jobParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_COMMAND, "command to submit counting job", true ));
    jobParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_COMMAND, "command to submit merging job", true ));
    jobParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_FILENAME, "filename to the couting job", true ));
    jobParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_FILENAME, "filename to the merging job", true ));
    jobParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_COUNT, "maximum number of simultaneous counting jobs", true ));
    jobParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_MERGE, "maximum number of simultaneous merging jobs", true ));

	getParser()->push_back(jobParser);
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
    Parameter (IProperties* props) : _props(props) {}
    //Simka& _simka;
    IProperties* _props;
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

	SimkaFusion<span> simkaAlgorithm (p._props);
	simkaAlgorithm.execute();

	/*
#ifdef SIMKA_MIN
	simkaAlgorithm.executeSimkamin();
#else
#endif*/
}};

void SimkaPotala::execute ()
{
	IProperties* input = getInput();
	//Parameter params(*this, getInput());
	Parameter params(input);

	size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    Integer::apply<Functor,Parameter> (kmerSize, params);
}





int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        SimkaPotala().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
