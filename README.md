
<details>
  <summary markdown="span">Table of contents</summary>

  [[_TOC_]]
</details>

# Simka & SimkaMin       

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


> This directory stores Simka and SimkaMin software. This readme focuses on Simka features. All information about SimkaMin is located in the [simkaMin](simkaMin/) directory. 


# Continuous integration status (master branch)

### Build status
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7-gitlab/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-debian7-64bits-gcc-4.7-gitlab/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1-gitlab/badge/icon)](https://ci.inria.fr/gatb-core/view/Simka/job/tool-simka-build-macos-10.9.5-gcc-4.2.1-gitlab/)


### SonarQube metrics

<details>
  <summary markdown="span">Click me to expand</summary>

  [![Lines of code](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=ncloc)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=ncloc)
  [![Comment line density](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=comment_lines_density)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=comment_lines_density)
  [![Coverage](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=coverage)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=coverage)

  [![Bugs](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=bugs)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=bugs)
  [![Vulnerabilities](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=vulnerabilities)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=vulnerabilities)
  [![Code Smells](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=code_smells)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=code_smells)

  [![New Bugs](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=new_bugs)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=new_bugs)
  [![New Vulnerabilities](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=new_vulnerabilities)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=new_vulnerabilities)
  [![New Code Smells](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=genscale:gatb:tools:simka:gitlab:master&metric=new_code_smells)](https://sonarqube.inria.fr/sonarqube/component_measures?id=genscale%3Agatb%3Atools%3Asimka%3Agitlab%3Amaster&metric=new_code_smells)
</details>

# What is Simka?

Simka is a de novo comparative metagenomics tool. Simka represents each dataset as a k-mer spectrum and compute several classical ecological distances between them.

Developper: [Gaëtan Benoit](http://people.rennes.inria.fr/Gaetan.Benoit/), PhD, former member of the [Genscale](http://team.inria.fr/genscale/) team at Inria.

Contact: claire dot lemaitre at inria dot fr

# References

* Simka: 
  Benoit G, Peterlongo P, Mariadassou M, Drezen E, Schbath S, Lavenier D, Lemaitre C. (2016) [Multiple comparative metagenomics using multiset k-mer counting](https://doi.org/10.7717/peerj-cs.94). PeerJ Computer Science 2:e94 

* SimkaMin:

  Gaetan Benoit, Mahendra Mariadassou, Stéphane Robin, Sophie Schbath, Pierre Peterlongo, Claire Lemaitre [SimkaMin: fast and resource frugal de novo comparative metagenomics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz685/5559271)
  Bioinformatics, https://doi.org/10.1093/bioinformatics/btz685

* Benoit G (2017) [Large scale de novo comparative metagenomics (PhD thesis in french)](https://tel.archives-ouvertes.fr/tel-01659395v2/).

# Install a binary release of simka

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

# Install simka from source code: git clone

Requirements: cmake 2.6+ and gcc 4.4.7+ (Linux) or clang 4.1+ (Mac OSX).

From the command-line:

    git clone https://github.com/GATB/simka.git
    cd simka
    sh INSTALL

See the INSTALL file for more information.

Then, you can try the software on your computer, as follows:

    cd example
    ./simple_test.sh

The installation creates 4 executables (./build/bin directory):

    simka: main software to be used for your analysis
    simkaCount, simkaMerge and simkaCountProcess: not to be used directly, called by 'simka'

All softwares must stay in the same folder; so, if you want to move them elsewhere on your system, consider to let them altogether.

For further instructions on using simka, see User Manual, below.

# Install simka from source code: using a source release archive

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

# Changelog

* version 1.5.1 Sept 05, 2019:   
    - simkaMin: easier usage of simkaMin, usefull for conda packaging 
* version 1.5 Jun 07, 2019:   
    - simkaMin software: faster results by subsampling the kmer space
* version 1.4 Jun 21, 2017:
	- update gatb-core to version 1.2.2
	- simka now provide gz compressed results
	- new scripts for result visualization
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

## Description
Simka computes several classical ecological distances between N (metagenomic) read sets based on k-mer counts.
Simka is implemented with the GATB library (http://gatb.inria.fr/).

## Input

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

## Output

### Temporary output

The option -out-tmp controls where the temporary files of Simka will be stored.

This option is mandatory since the disk usage of Simka can be high depending on the input size.

This option must target a directory on your faster disk with some free space.

One may want to add new datasets to existing Simka results without recomputing everything again (for instance, if your metagenomic project is incomplete).
This can only be achieved by keeping those temporary files on the disk using the option -keep-tmp of Simka.

### Result output

Simka results are distance matrices. A distance matrix is a squared matrix of size N (where N is the number of input datasets). Each value in the matrix give you the distance between a pair of datasets. These values are usually in the range [0, 1]. A distance value of 0 means that the pair of dataset is perfectly similar. The higher the distance value is, the more dissimilar is the pair of datasets.

Simka results will be stored in the directory indicated by -out option.

By default, Simka compute an abundance-based Bray-Curtis distance matrix and a presence-absence-based Jaccard distance matrix.

The option -simple-dist allows to compute more ecology distances which are fast to compute (Chord, Hellinger, Kulczinski...).

The option -complex-dist allows to compute others ecology distances which can be very long to compute (Jensen-Shannon, Canberra, Whittaker...).

The matrice names follow this template:

    mat_[abundance|presenceAbsence]_[distanceName].csv.gz

The distance matrices containing ‘simka’ are distances introduces by the comparead method.
These distances have the advantage of having a symmetrical and asymmetrical version.

## Visualize simka results

Simka results can be visualized through heatmaps, hierarchical clustering and PCA (MDS or PCoA to be exact).
	
Requirements: R, gplots package (only for heatmap)

Use the script run-visualization.py (located in "scripts/visualization" folder).

Example: 

```bash
python run-visualization.py -in simka_results_dir -out output_figures_dir -pca -heatmap -tree
```

where simka_results_dir is the folder containing the distances matrices of Simka (-out)

Figures can be annotated by providing a metadata data in standard csv format:

```bash
DATASET_ID;VARIABLE_NAME_1;VARIABLE_NAME_2
A;1;aquatic
B;1;human
C;2;human
D;2;soil
E;3;soil
```

An example of this table is given at ./example/dataset_metadata.csv

Dataset ID in the metadata table must match with the dataset ID in simka distance matrices

Add the following options to activate annotations:

```bash
-metadata-in: filename to a metadata table
-metadata-variable: the name of the variable that you want to display in figures (the name of the column), for instance VARIABLE_NAME_1 in example above
```

Visualization example commands are given when running simka example (./example/simple_test.sh).

## Usage for simka

To see simka in-line help:

```bash
./bin/simka
```


## Simka command examples

Run the toy example:

```bash
./bin/simka -in example/simka_input.txt -out results -out-tmp temp_output
```

Compute all the distances that Simka can provide (Bray-Curtis, Jensen-Shannon…):

```bash
./bin/simka … -simple-dist -complex-dist
```

Change the kmer size

```bash
./bin/simka … -kmer-size 31
```

Filter kmers seen one time (potentially erroneous) and very high abundance kmers (potentially contaminants):

```bash
./bin/simka … -abundance-min 2 -abundance-max 200
```

Filter over the sequences of the reads and k-mers:

Minimum read size of 90. Discards low complexity reads and k-mers (shannon index < 1.5)

```bash
./bin/simka … -min-read-size 90 -read-shannon-index 1.5 -kmer-shannon-index 1.5
```

Consider a subset of the reads of the input dataset (for dataset with non-uniform reads per sample):

Considers all the reads of each samples (default)

```bash
./bin/simka … -max-reads -1
```

Let Simka compute automatically the maximum of read per samples (normalization)

```bash
./bin/simka … -max-reads 0
```

Used only the first 1000 reads of each samples:

```bash
./bin/simka … -max-reads 1000
```

Allow more memory and cores improve the execution time:

```bash
./bin/simka … -max-memory 20000 -nb-cores 8
```


## Computer cluster options

Simka can be ran on computer cluster equipped of a job scheduling system such as SGE. Giving a job file template and a submission command, Simka will take care of creating and synchronizing the jobs until the end of the execution.

You must provide the filenames to two job templates, one for counting and one for merging (-count-file -count-merge).

There are example of file templates in the folder ‘example/potara_job’.

And you must provide a submission command for both job (-count-cmd -merge-cmd)

Example for SGE:

```bash
-count-cmd ‘qsub  -pe make 8’ -merge-cmd qsub
```

The option -max-count and -max-merge controls the maximum of simultaneous jobs. They have to be fixed if you system have a maximum of jobs restriction.

Command example:

```bash
./bin/simka … -count-file example/potara_job/sge/job_count -merge-file example/potara_job/sge/job_merge \
-count-cmd qsub -pe make 34 -merge-cmd qsub \
-max-count 6 -max-merge 18 -nb-cores 200 -max-memory 500000
```

Simka will run a maximum of 6 simultaneous counting jobs, each using 200/6 cores and 500000/6 MB of memory. Simka will run a maximum of 18 merging jobs. A merging job can not be ran on more than 1 core and use very low memory. By default Simka use -nb-cores/2 counting jobs simultaneously and -nb-cores merging jobs simultaneously.


## Possible issues with Simka

### TOO MUCH OPENED FILES

Simka is a disk-based method. Depending on the chosen options (-nb-cores -max-memory), it is possible that Simka required a lot of open files.

You can fix this issue in two ways:

* increasing the maximum open files limit imposed by your system: ulimit -n maxFiles
* reducing the number of files opened by Simka by using the option -max-count and -max-merge
