--------------------------------------------------------------------------------
# RELEASE 1.1.1

* Re-design to support variable number of kmer sizes
 => now, one can use the cmake variable KSIZE_LIST, for instance "cmake -DKSIZE_LIST="32 64 96" ..

* Allows "auto" value for the -abundance-min parameter

--------------------------------------------------------------------------------
# RELEASE 1.1.0

* Re-design of the SortingCountAlgorithm with introduction of interface ICountProcessor
 => it should allow development of new tools based on kmers counting

--------------------------------------------------------------------------------
# RELEASE 1.0.8

* Correction of memory alignment issue on MacOs in some cases

* Re-introduce multi-passes management in DSK

* Correction of passes number configuration with some banks inputs

* Temporary files have now unique names so dbgh5 can be launched several times in the same working directory

--------------------------------------------------------------------------------
# RELEASE 1.0.7

* Correction of scripts for new project creation and delivery process

--------------------------------------------------------------------------------
# RELEASE 1.0.6

* Speed up from x2 to x3 for kmer counting and graph construction phases (optimizations based on minimizers and improved Bloom filters). GATB's k-mer counter has been improved using techniques from KMC2, to achieve competitive running times compared to KMC2.

* Ability to store arbitrary information associated to each kmer of the graph, enabled by a minimal perfect hash function (costs only 2.61 bits/kmer of memory)

* Improved API with new possibilities (banks and kmers management)

* Many new snippets showing how to use the library.


--------------------------------------------------------------------------------
# RELEASE 1.0.5

Modifications of Kmer::Model class for kmers management
* better implementation (factorization, optimization)
* introduction of minimizers concept

WARNING ! These modifications introduced small API changes. Please read snippets kmer2 and kmer5 to see how to handle kmers now.

