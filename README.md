# Simka
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is Simka?

Simka is a comparative metagenomics method dedicated to NGS datasets. It computes a large collection of distances classically used in ecology to compare communities by approximating species counts by k-mer counts.

Developper: Gaëtan Benoit.
Contact: gaetan.benoit@inria.fr

#Reference
	
G. Benoit, P. Peterlongo, D. Lavenier, C. Lemaitre. (2015) [Simka: fast kmer-based method for estimating the similarity between numerous metagenomic datasets](https://hal.inria.fr/hal-01180603). Hal-Inria

G. Benoit, P. Peterlongo, M. Mariadassou, E Drezen, S. Schbath, D. Lavenier, C. Lemaitre. Multiple comparative metagenomics using multiset k-mer counting. *Submitted*.

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

    cd example
    ./simple_test.sh

The installation creates 3 executables (./build/bin directory):

    simka: main software to be used for your analysis
    simkaCount: not to be used directly, called by 'simka'
    simkaMerge: not to be used directly, called by 'simka'

All softwares must stay in the same folder; so, if you want to move them elsewhere on your system, consider to let them altogether.

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

    cd example
    ./simple_test.sh

For further instructions on using simka, see User Manual, below.

#Changelog

* version 1.0.1  March 16, 2016: minor updates ang bug fixes, first release on Github
* version 1  Feb 16, 2016: stable version
* version 0.1  May 28, 2015: initial public release

# User manual

##Description
Simka computes several ecology distances between N metagenomic read sets at the k-mer level.
Simka is implemented with the GATB library (http://gatb.inria.fr/).


##Input
	
The input file (-in) lists the samples. This sample can be in fasta, fastq and in gzip compressed format (.gz).

One sample per line with the following syntax:

    ID1: filename.fasta
    ID2: filename.fasta
    ID3: filename.fasta

The sample ID in the name that will appear in the headers of the distance matrices.

Simka supports the management of paired and concatenated files.

If you have paired files, you can list them separated by a ‘;’:

    ID1_paired: filename_pair1.fasta ; filename_pair2.fasta

If you want to concatenate some files (for example, if you have multiple experiment for a single environment), you can list them separated by a ‘,’:

    ID1_concat: filenameX.fasta.gz , filenameY.fasta.gz

You can combine concatenated and paired operations:

    ID1_concat: filenameX_1.fasta.gz , filenameY_1.fasta.gz ; filenameX_2.fasta.gz , filenameY_2.fasta.gz

There is a difference between paired and concatenated files depending on the value of the option -max-reads.

Example:

If -max-reads is set to 100, then Simka will considered the 100 first reads of the first paired files and the 100 first reads of the second paired files…

##Output

###Temporary output

The option -out-tmp controls where the temporary files of Simka will be stored.

This option is mandatory since the disk usage of Simka can be high depending on the input size.

This option must target a directory on your faster disk with some free space.

At the end of an execution, Simka does not remove the temporary files.
These files can be re-used in case you want to add new samples in the input files.
In this case, Simka will not count again the samples already counted.

###Results output

The option -out is the directory which will hold the distances matrices.

By default, Simka provides a abundance-based Bray-Curtis distance (and Jaccard) and all presence-absence-based distances.

The option -simple-dist allow to compute more ecology distances which are fast to compute (Chord, Hellinger, Kulczinski).

The option -complex-dist allow to compute others ecology distances which are long to compute (Jensen-Shannon, Canberra, Whittaker).

The matrice names follow this template:

    mat_[abundance|presenceAbsence]_[distanceName].csv

The distance matrices containing ‘simka’ are distances introduces by the comparead method (See equation 1 of the Simka paper).
These distances have the advantage of having a symmetrical and asymmetrical version.

##Generating heatmaps and clustering
	
Requirements: python,R and gplots package

Run the script create_heatmaps.py (located in "scripts" folder) in the scripts folder.

Example: 

    python create_heatmaps.py matricesFolder

where matricesFolder in the folder containing the distances matrices of Simka (-out)


##Usage for simka

To see simka in-line help:

    ./bin/simka


##The options of simka by command examples

Run the toy example:

    ./bin/simka -in example/simka_input.txt -out results -out-tmp temp_output

Compute all the distances that Simka can provide (Bray-Curtis, Jensen-Shannon…):

    ./bin/simka … -simple-dist -complex-dist

Change the kmer size

    ./bin/simka … -kmer-size 31

Filter kmers seen one time (potentially erroneous) and very high abundance kmers (potentially contaminants):

    ./bin/simka … -abundance-min 2 -abundance-max 200

Filter over the sequences of the reads and k-mers:

Minimum read size of 90. Discards low complexity reads and k-mers (shannon index < 1.5)

    ./bin/simka … -min-read-size 90 -read-shannon-index 1.5 -kmer-shannon-index 1.5

Consider a subset of the reads of the input dataset (for dataset with non-uniform reads per sample):

Considers all the reads of each samples (default)

    ./bin/simka … -max-reads -1

Let Simka compute automatically the maximum of read per samples (normalization)

    ./bin/simka … -max-reads 0

Used only the first 1000 reads of each samples:

    ./bin/simka … -max-reads 1000

Allow more memory and cores improve the execution time:

    ./bin/simka … -max-memory 20000 -nb-cores 8


##Computer cluster options

Simka can be ran on computer cluster equipped of a job scheduling system such as SGE. Giving a job file template and a submission command, Simka will take care of creating and synchronizing the jobs until the end of the execution.

You must provide the filenames to two job templates, one for counting and one for merging (-count-file -count-merge).

There are example of file templates in the folder ‘example/potara_job’.

And you must provide a submission command for both job (-count-cmd -merge-cmd)

Example for SGE:

    -count-cmd ‘qsub  -pe make 8’ -merge-cmd qsub

The option -max-count and -max-merge controls the maximum of simultaneous jobs. They have to be fixed if you system have a maximum of jobs restriction.

Command example:

    ./bin/simka … -count-file example/potara_job/sge/job_count -merge-file example/potara_job/sge/job_merge \
    -count-cmd qsub -pe make 34 -merge-cmd qsub \
    -max-count 6 -max-merge 18 -nb-cores 200 -max-memory 500000
    
Simka will run a maximum of 6 simultaneous counting jobs, each using 200/6 cores and 500000/6 MB of memory. Simka will run a maximum of 18 merging jobs. A merging job can not be ran on more than 1 core and use very low memory. By default Simka use -nb-cores/2 counting jobs simultaneously and -nb-cores merging jobs simultaneously.


##Possible issues with Simka

TOO MUCH OPENED FILES

Simka is a disk-based method. Depending on the chosen options (-nb-cores -max-memory), it is possible that Simka required a lot of open files.

You can fix this issue in two ways:

* increasing the maximum open files limit imposed by your system: ulimit -n maxFiles
* reducing the number of files opened by Simka by using the option -max-count and -max-merge
