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

#include "Simka.hpp"
#include "SimkaAlgorithm.hpp"


IOptionsParser* Simka::createOptionsParser (IOptionsParser* parent)
{
    IOptionsParser* parser = parent; //new OptionsParser ("Simka");

	//Main parser
    parser->push_front (new OptionNoParam (STR_SIMKA_COMPUTE_DATA_INFO, "compute (and display) information before running Simka, such as the number of reads per dataset", false));
    parser->push_front (new OptionNoParam (STR_SIMKA_KEEP_TMP_FILES, "keep temporary files", false));
    parser->push_front (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
    parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "output directory for result files (distance matrices)", false, "./simka_results"));
    parser->push_front (new OptionOneParam (STR_URI_INPUT, "input file of samples. One sample per line: id1: filename1...", true));


    //parser->push_back (new OptionOneParam (STR_URI_OUTPUT_TMP, "output directory for temporary files", true));
	//IOptionsParser* parser = getParser();
	//IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();
	//parser->push_back(dskParser);
	//dskParser->setVisible(false);
	//cout << parser->getParser(STR_NB_CORES) << endl;
	parser->getParser(STR_NB_CORES)->setVisible(false);

	//parser->push_back(new OptionOneParam(parser->getParser(STR_NB_CORES)->getName(), parser->getParser(STR_NB_CORES)->getHelp(), false, "0"));
	//parser->push_front(dskParser->getParser (STR_URI_OUTPUT_TMP));
	//dskParser->getParser (STR_URI_OUTPUT_TMP)->setMandatory
    //parser->push_front(dskParser->getParser (STR_URI_OUTPUT));
    //parser->getParser (STR_URI_OUTPUT)->setHelp("output directory for result files (similarity matrix, heatmaps)");
    //parser->push_front(dskParser->getParser (STR_URI_INPUT));
    //parser->getParser(STR_URI_INPUT)->setHelp("input file of datasets. One dataset per line: id filename1 filename2...");

    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_URI_OUTPUT_TMP)))  {  p->s; }

	//Distance parser
    IOptionsParser* distanceParser = new OptionsParser ("distance");
    distanceParser->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES, "compute all simple distances (Chord, Hellinger...)", false));
    distanceParser->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES, "compute all complex distances (Jensen-Shannon...)", false));


	//Kmer parser
    IOptionsParser* kmerParser = new OptionsParser ("kmer");
    kmerParser->push_back (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "21"));
    //kmerParser->push_back(dskParser->getParser (STR_KMER_SIZE));
    //kmerParser->push_back(new OptionOneParam (STR_KMER_PER_READ.c_str(), "number of selected kmers per read", false, "0"));
    //kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "1"));
    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "min abundance a kmer need to be considered", false, "2"));
    kmerParser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "max abundance a kmer can have to be considered", false, "999999999"));
    //kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MIN));
    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    //if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

    //kmerParser->push_back(dskParser->getParser (STR_KMER_ABUNDANCE_MAX));
    //kmerParser->push_back(dskParser->getParser (STR_SOLIDITY_KIND));
    //kmerParser->getParser (STR_SOLIDITY_KIND)->setHelp("TODO");
    //kmerParser->push_back (new OptionNoParam (STR_SIMKA_SOLIDITY_PER_DATASET.c_str(), "do not take into consideration multi-counting when determining solid kmers", false ));
    kmerParser->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX.c_str(), "minimal Shannon index a kmer should have to be kept. Float in [0,2]", false, "0" ));


    //Read filter parser
    IOptionsParser* readParser = new OptionsParser ("read");
    readParser->push_back (new OptionOneParam (STR_SIMKA_MAX_READS.c_str(), "maximum number of reads per sample to process. Can be -1: use all reads. Can be 0: estimate it", false, "-1" ));
    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE.c_str(), "minimal size a read should have to be kept", false, "0" ));
    readParser->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX.c_str(), "minimal Shannon index a read should have to be kept. Float in [0,2]", false, "0" ));

    //Core parser
    IOptionsParser* coreParser = new OptionsParser ("core");
    coreParser->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "0"));
    coreParser->push_back (new OptionOneParam (STR_MAX_MEMORY, "max memory (MB)", false, "5000"));
    //coreParser->push_back(dskParser->getParser ());
    //coreParser->push_back(dskParser->getParser (STR_MAX_DISK));

    //Distances
    //IOptionsParser* distanceParser = new OptionsParser ("distances");
    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_BRAYCURTIS.c_str(), "compute Bray Curtis distance"));
    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_CHORD.c_str(), "compute Chord distance"));
    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_HELLINGER.c_str(), "compute Hellinger distance"));
    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_CANBERRA.c_str(), "compute Canberra distance"));
    //distanceParser->push_back (new OptionNoParam (STR_SIMKA_DISTANCE_KULCZYNSKI.c_str(), "compute Kulczynski distance"));


	parser->push_back(distanceParser);
	parser->push_back(kmerParser);
	parser->push_back(readParser);
	parser->push_back(coreParser);
	//parser->push_back(distanceParser);


	IOptionsParser* dskParser = SortingCountAlgorithm<>::getOptionsParser();

    if (Option* p = dynamic_cast<Option*> (dskParser->getParser(STR_MINIMIZER_SIZE)))  {  p->setDefaultValue ("7"); }
	parser->push_back(dskParser);
	dskParser->setVisible(false);
    if (Option* p = dynamic_cast<Option*> (parser->getParser(STR_SOLIDITY_KIND)))  {  p->setDefaultValue ("all"); }

    return parser;
}


Simka::Simka()  : Tool ("Simka")
{

	Simka::createOptionsParser(getParser());

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
	SimkaAlgorithm<span> simkaAlgorithm (p._props);
	simkaAlgorithm.execute();

	/*
#ifdef SIMKA_MIN
	simkaAlgorithm.executeSimkamin();
#else
#endif*/
}};

void Simka::execute ()
{
	IProperties* input = getInput();

	//Parameter params(*this, getInput());
	Parameter params(input);

	size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);
	/*
	params._kmerSize = getInput()->getInt (STR_KMER_SIZE);
	params._inputFilename = input->getStr(STR_URI_INPUT);
	params._outputDir = input->get(STR_URI_OUTPUT) ? input->getStr(STR_URI_OUTPUT) : "./";
	params._abundanceThreshold.first = input->getInt(STR_KMER_ABUNDANCE_MIN);
	params._abundanceThreshold.second = input->getInt(STR_KMER_ABUNDANCE_MAX);
	params._soliditySingle = input->get(Simka::STR_SOLIDITY_PER_DATASET);

	cout << params._soliditySingle << endl;
*/


    /** We launch the tool with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<Functor,Parameter> (kmerSize, params);
}



