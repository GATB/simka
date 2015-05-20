

#include "Simka.hpp"
#include "SimkaAlgorithm.hpp"

Simka::Simka()  : Tool ("Simka")
{

	IOptionsParser* parser = getParser();

	IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();

    //if (IOptionsParser* input = dskParser->getParser (STR_KMER_ABUNDANCE_MIN_THRESHOLD))  {  input->setVisible (false);  }

	//dskParser->
	parser->push_back (dskParser, 1);

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
    Parameter (Simka& simka, IProperties* props) : simka(simka), props(props) {}
    Simka&       simka;
    IProperties* props;
};

template<size_t span> struct Functor  {  void operator ()  (Parameter parameter)
{
	SimkaAlgorithm<span> simkaAlgorithm (parameter.props);
	simkaAlgorithm.execute();
}};

void Simka::execute ()
{
    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** We launch the tool with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<Functor,Parameter> (kmerSize, Parameter (*this, getInput()));
}



