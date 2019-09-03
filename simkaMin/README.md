# SimkaMin
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## What is SimkaMin?

As in the case of Simka, SimkaMin is a *de novo* comparative metagenomics tool. The difference with Simka stands in the fact that SimkaMin outputs approximate (but very similar) results by subsampling the kmer space. With this strategy, and with default parameters, SimkaMin is an order of magnitude faster, uses 10 times less memory and 70 times less disk than Simka. 

Developper: [Gaëtan Benoit](http://people.rennes.inria.fr/Gaetan.Benoit/), PhD, former member of the [Genscale](http://team.inria.fr/genscale/) team at Inria.

Contact: claire dot lemaitre at inria dot fr

## References

Benoit G,  Mariadassou M, Robin S, Schbath S, Peterlongo P and Lemaitre C. (2019) [SimkaMin: fast and resource frugal *de novo* comparative metagenomics](https://doi.org/10.1093/bioinformatics/btz685). Bioinformatics

Benoit G, Peterlongo P, Mariadassou M, Drezen E, Schbath S, Lavenier D, Lemaitre C. (2016) [Multiple comparative metagenomics using multiset k-mer counting](https://doi.org/10.7717/peerj-cs.94). PeerJ Computer Science 2:e94 

Benoit G (2017) [Large scale de novo comparative metagenomics (PhD thesis in french)](https://tel.archives-ouvertes.fr/tel-01659395v2/).

## Install simkaMin

SimkaMin comes with Simka installation. Refer to [Simka install instructions](../README.md). 

## User manual

### Description
SimkaMin computes Bray-Curtis (abundance based) and Jaccard (presence/absence based) distances between N (metagenomic) read sets based on subsamples of k-mer counts.

Basically it takes as input the N metagenomic read sets and it outputs two matrices respectively providing the pairwise Bray-Curtis and the Jaccard distances between each dataset pairs. 

### A simple command example

Run the toy example:

```bash
./simkaMin/simkaMin.py -in example/simka_input.txt -out results 
```


### Input

The input file (`-in`) lists the datasets. These datasets can be in fasta, fastq and in gzip compressed format (.gz).

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

Paired syntax is only usefull if the `-max-reads` option of SimkaMin is set.

Example:

If `-max-reads` is set to 100, then Simka will considered the 100 first reads of the first paired files and the 100 first reads of the second paired files…

### Output

SimkaMin results are an abundance-based Bray-Curtis distance matrix `mat_presenceAbsence_jaccard.csv.gz` and a presence-absence-based Jaccard distance matrix `mat_abundance_braycurtis.csv.gz`. A distance matrix is a squared matrix of size N (where N is the number of input datasets). Each value in the matrix gives the distance between a pair of datasets. These values are in the range [0, 1]. A distance value of 0 means that the pair of datasets is perfectly similar. The greater the distance value is, the more dissimilar is the pair of datasets.

SimkaMin results will be stored in the directory indicated by `-out` option.

#### Visualize SimkaMin results

SimkaMin results can be visualised through heatmaps, hierarchical clustering and PCA. This module is common with the Simka visualisation script `../scripts/visualization/run-visualization.py`.

Please refer to the documentation provided in the [Simka Readme file](../README.md). 	


## Usage

To see simka in-line help:

```bash
./simkaMin/simkaMin.py 
```


## Simka command examples

Run the toy example:

```bash
./simkaMin/simkaMin.py -in example/simka_input.txt -out results 
```

Change the kmer size

```bash
./simkaMin/simkaMin.py … -kmer-size 31
```

Change the sub-sampling effort (default 1 million kmers are used per read set)

```bash
./simkaMin/simkaMin.py … -nb-kmers 10000
```

Filter kmers seen one time (potentially erroneous):

```bash
./simkaMin/simkaMin.py … -filter
```

Consider all the reads of each samples (set 0 to use all reads)

```bash
./simkaMin/simkaMin.py … -max-reads 0
```

Use only the first 1000 reads of each sample:

```bash
./simkaMin/simkaMin.py … -max-reads 1000
```

Allow more memory and cores to improve the execution time:

```bash
./simkaMin/simkaMin.py … -max-memory 20000 -nb-cores 8
```

Filter low complexity reads

```bash
./simkaMin/simkaMin.py … -min-shannon-index 1
```

Filter small reads 

```bash
./simkaMin/simkaMin.py … -min-read-size 80
```

Update existing results with additional datasets

```bash
./simkaMin/simkaMin_update.py -in another_simka_input.txt -in-to-update results/simkamin
# updated matrices will be in dir results/simkamin/
```



