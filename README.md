# Simka
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is Simka?

Simka is a de novo comparative metagenomics tool. Simka represents each dataset as a k-mer spectrum and compute several classical ecological distances between them.

Developper: Gaëtan Benoit.

Contact: gaetan.benoit@inria.fr

# Reference
	
Benoit G, Peterlongo P, Mariadassou M, Drezen E, Schbath S, Lavenier D, Lemaitre C. (2016) [Multiple comparative metagenomics using multiset k-mer counting](https://doi.org/10.7717/peerj-cs.94). PeerJ Computer Science 2:e94 

Benoit G, Peterlongo P, Lavenier D, Lemaitre C. (2015) [Simka: fast kmer-based method for estimating the similarity between numerous metagenomic datasets](https://hal.inria.fr/hal-01180603). Hal-Inria

#Install a binary release of simka

Retrieve the binary archive file from one of the official simka releases (see "Releases" tab on the Github web page of the simka project); file name is "simka-xyz-bin-Darwin.tar.gz" or "simka-xyz-bin-Linux.tar.gz" (where xyz is a release version).

Then, from the command-line:

    gunzip simka-xyz-bin-Dawrin.tar.gz
    tar -xf simka-xyz-bin-Dawrin.tar
    cd simka-xyz-Dawrin
    chmod +x bin/* example/*.sh

Binary of simka is in folder "bin". You can try the software on your computer, as follows:

    cd example
    ./simple_test.sh

In case the software does not run appropriately on your system, you should consider to install it from its source code, as explained below.

For further instructions on using simka, see User Manual, below.

#Install simka from source code: git clone

Requirements: cmake 2.6+ and gcc 4.5+ (Linux) or clang 4.1+ (Mac OSX).

From the command-line:

    git clone https://github.com/GATB/simka.git
    cd simka
    sh INSTALL

See the INSTALL file for more information.

Then, you can try the software on your computer, as follows:

    python ./example/1-basic_usage/1-simple_test.py

For further instructions on using simka, see User Manual, below.

#Install simka from source code: using a source release archive

Requirements: cmake 2.6+ and gcc 4.5+ (Linux) or clang 4.1+ (Mac OSX).

Retrieve the source code archive file from one of the official simka releases (see "Releases" tab on the Github web page of the simka project); file name is "simka-xyz-Source.tar.gz" (where xyz is a release version).

Then, from the command-line:

    gunzip simka-xyz-Source.tar.gz
    tar -xf simka-xyz-Source.tar
    cd simka-xyz-Source
    sh INSTALL

Then, you can try the software on your computer, as follows:

    python ./example/1-basic_usage/1-simple_test.py

For further instructions on using simka, see User Manual, below.


#Changelog

* version 2.0.0 Nov 28, 2016:
	- simka code has been refactored for robustness and flexibility
	- existing run of simka can be updated without recomputing everything
	- simka now provides compressed results (.gz)
* version 1.3.2 Oct 25, 2016:
	- improve memory usage of symetrical distances
	- option -data-info to compute information on the input data (nb reads per dataset...)
	- intermediate merge sort passes to handle large number of datasets
	- prevent distances from producing nan value
	- fix bug that occur during k-mer counting
* version 1.3.0 July 29, 2016:
	- Bray-Crutis computed by default
	- Better k-mer statistics
	- Fix bug in script for creating heatmaps
	- Add "all in memory" k-mer counter when k <= 15
	- Fine grain paralellization for computing distances
	- Clean all memory leaks with valgrind
	- Update help messages
	- Redirect stdout and stderr of parallel processes in specific log files 
* version 1.0.1  March 16, 2016: minor updates ang bug fixes, first release on Github
* version 1  Feb 16, 2016: stable version
* version 0.1  May 28, 2015: initial public release

# User manual

##Description
Simka computes several ecology distances between N (metagenomic) read sets at the k-mer level.
Simka is implemented with the GATB library (http://gatb.inria.fr/).

##Simka overall organisation
* example
	- 1-basic_usage (learn how to use simka by running simple examples)
	- 2-hpc_usage (learn how to run simka on cloud or grid systems)
	- 3-simka-pipeline (understand how works each piece of Simka)
	- data (show how to layout simka input)
* scripts
	- visualization (collection of R scripts to visualize Simka results)
	- simka2
		- bin (simka binaries location after compilation)
		- core (collection of python scripts for running Simka)
		- simka.py (main script for running simka)
		- simka-hpc.py (script for running simka in HPC mode)
* src
	- source code of Simka written in c++
* tests
	- validation tests of Simka

##Input
	
The input file (-in) lists the datasets. These datasets can be in fasta, fastq and in gzip compressed format (.gz).

One dataset per line with the following syntax (you can put any number of spaces and/or tabs between syntax):

    ID1: filename.fasta
    ID2: filename.fasta
    ID3: filename.fasta

The dataset ID in the name that will appear in the headers of the distance matrices.

You can find a simka input file in example directory: ./example/data/simka_input.txt

If a given datset has been splitted in several parts, Simka can automatically concatenate them.

    ID1: filename_part1.fasta , filename_part2.fasta , ...
    
If you have paired files, you can list them separated by a ‘;’:

    ID1: filename_pair1.fasta ; filename_pair2.fasta

You can combine concatenated and paired operations:

    ID1: filename_part1_pair1.fasta , filename_part2_pair1.fasta ; filename_part1_pair2.fasta , filename_part2_pair2.fasta

Paired syntax is only usefull if the -max-reads option of Simka is set.

Example:

If -max-reads is set to 100, then Simka will considered the 100 first reads of the first paired files and the 100 first reads of the second paired files…

##Output

###Temporary output

The option -out-tmp controls where the temporary files of Simka will be stored.

This option is mandatory since the disk usage of Simka can be high depending on the input size.

This option must target a directory on your faster disk with some free space.

One may want to add new datasets to existing Simka results without recomputing everything again (for instance, if your metagenomic project is incomplete).
This can only be achieved by keeping those temporary files on the disk using the option -keep-tmp of Simka.

###Result output

Simka results are distance matrices. A distance matrix is a squared matrix of size N (where N is the number of input datasets). Each value in the matrix give you the distance between a pair of datasets. These values are usually in the range [0, 1]. A distance value of 0 means that the pair of dataset is perfectly similar. The higher the distance value is, the more dissimilar is the pair of datasets.

Simka results will be stored in the directory indicated by -out option.

By default, Simka compute a abundance-based Bray-Curtis distance matrix and a presence-absence-based Jaccard distance matrix.

The option -simple-dist allows to compute more ecology distances which are fast to compute (Chord, Hellinger, Kulczinski...).

The option -complex-dist allows to compute others ecology distances which can be very long to compute (Jensen-Shannon, Canberra, Whittaker...).

The matrice names follow this template:

    mat_[abundance|presenceAbsence]_[distanceName].csv.gz

The distance matrices containing ‘simka’ are distances introduces by the comparead method.
These distances have the advantage of having a symmetrical and asymmetrical version.

The result directory also contains a directory named "matrix_binary" that can be used by the simka distance exporter to create quickly new distance matrices from supplied list of dataset ID (see "example" section).

##Visualize simka results

Simka results can be visualized through heatmaps, hierarchical clustering and PCA (MDS or PCoA to be exact).
	
Requirements: R, gplots package (only for heatmap)

Use the script run-visualization.py (located in "scripts/visualization" folder).

Example: 

    python run-visualization.py -in simka_results_dir -out output_figures_dir -pca -heatmap -tree

where simka_results_dir is the folder containing the distances matrices of Simka (-out)

To learn advanced usage, run example: ./example/1-basic_usage/2-visualization.py.
For instance, you can add annotations to figures (colors) by provding a metadata table in a specific format.

##Usage for simka

To see simka in-line help:

    python ./scripts/simka2/simka.py

##Simka command examples

Commands listed below are simple and basic usage of Simka.
To learn quickly and easily how to use simka and discover advanced usage, run python scripts in "example" directory.

Run the toy example:

    python ./scripts/simka2/simka.py -in ./example/data/simka_input.txt -out simka_results -out-tmp simka_temp

Compute all the distances that Simka can provide (Bray-Curtis, Jensen-Shannon…):

    python ./scripts/simka2/simka.py … -simple-dist -complex-dist

Change the kmer size

    python ./scripts/simka2/simka.py … -kmer-size 31

Filter kmers seen one time (potentially erroneous) and very abundant kmers (potentially contaminants):

    python ./scripts/simka2/simka.py … -abundance-min 2 -abundance-max 10000

Control read quality (minimum read size of 90 and discards low complexity reads)

    python ./scripts/simka2/simka.py … -min-read-size 90 -read-shannon-index 1.5

Consider a subset of the reads of the input dataset (for dataset with non-uniform reads per sample):

Considers all the reads of each samples (default)

    python ./scripts/simka2/simka.py … -max-reads 0

Use only the first 1000 reads of each samples:

    python ./scripts/simka2/simka.py … -max-reads 1000

Allow more memory and cores improve the execution time:

    python ./scripts/simka2/simka.py … -max-memory 20000 -nb-cores 8

Run simka in HPC mode, with a maximum of 50 jobs simultaneously, 4 cores and 4 GB of memory per job:
(The job submit command is specific to your HPC system)

    python ./scripts/simka2/simka-hpc.py … -nb-cores 4 -max-memory 4000 -max-jobs 50 -submit-command "qsub -pe make 4 -M 4000"
    
##Simka examples

Simka examples are located in the "example" dir. This folder is organised as follow:
* 1-basic_usage
	- learn how to use simka by running simple examples
* 2-hpc_usage
	- learn how to run Simka with high performance computing (HPC)
* 3-simka-pipeline
	- understand how works each piece of Simka
* data
	- show how to layout simka input

A brief description of each example:
* 1-basic_usage
	- 1-simple_test.py
		- run Simka with its default parameters
	- 2-visualization.py
		- run a collection of R scripts for visualizing simka results
	- 3-export_distances.py
		- extract quickly columns and rows from the distance matrices
		- reorder dataset IDs in distance matrices
	- 4-add_new_datasets.py
		- show how it is possible to add new datasets to existing distance matrices
	- 5-simka_parameters.py
		- overview of simka optional parameters
* 2-hpc_usage
	- 1-submit_with_job_command.py
		- for system that submit job through a submit command
	- 2-submit_with_job_file.py
		- for system that ssubmit job through a command and a job file
* 3-simka_pipeline
	- 1-init_simka_database.py
		- database initilization (records simka parameters like k-mer size)
	- 2-compute_kmers_spectrums.py
		- compute k-mer spectrums (one per dataset)
	- 3-compute_distances.py
		- compute distances between k-mer spectrums


