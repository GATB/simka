

import os, sys, argparse
from core.simka2_utils import ArgumentFormatterSimka, SimkaParser


parser = SimkaParser(formatter_class=ArgumentFormatterSimka)


parserMain = parser.add_argument_group("main options")
parserCore = parser.add_argument_group("HPC options")
parserDistance = parser.add_argument_group("distance options")
parserKmer = parser.add_argument_group("k-mer options")
parserRead = parser.add_argument_group("read options")
parserDev = parser.add_argument_group("advanced (developer) options")

parserMain.add_argument('-in', action="store", dest="input_filename", help="input file of samples. One sample per line: id1: filename1...", required=True)
parserMain.add_argument('-out-tmp', action="store", dest="output_dir_temp", help="output directory for temporary files", required=True)
#parserMain.add_argument('-simka-bin', action="store", dest="simka_bin_dir", help="dir containing simka binaries", required=True)
parserMain.add_argument('-out', action="store", dest="output_dir", default="./simka_results", help="output directory for result files (distance matrices)")
parserMain.add_argument('-keep-tmp', action="store_true", dest="keep_tmp", help="keep temporary files. Allow to update existing run of simka")

parserDistance.add_argument('-simple-dist', action="store_true", dest="simple_dist", help="compute all simple distances (Chord, Hellinger...)")
parserDistance.add_argument('-complex-dist', action="store_true", dest="complex_dist", help="compute all complex distances (Jensen-Shannon...)")

parserKmer.add_argument('-kmer-size', action="store", dest="kmer_size", help="size of a kmer", default="21")
parserKmer.add_argument('-abundance-min', action="store", dest="abundance_min", help="min abundance a kmer need to be considered", default="2")
parserKmer.add_argument('-abundance-max', action="store", dest="abundance_max", help="max abundance a kmer can have to be considered", default="0")
#parser.add_argument('-kmer-shannon-index', action="store", dest="simple_dist")

parserRead.add_argument('-max-reads', action="store", dest="max_reads", default="0", help="maximum number of reads per sample to process")
parserRead.add_argument('-min-read-size', action="store", dest="min_read_size", default="0", help="minimal size a read should have to be kept")
parserRead.add_argument('-min-shannon-index', action="store", dest="min_read_shannon_index", default="0", help="minimal Shannon index a read should have to be kept. Float in [0,2]")

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores per job", required=True)
parserCore.add_argument('-max-memory', action="store", dest="max_memory", help="max memory (MB) per job", required=True)
parserCore.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of jobs that can be submitted simultaneously", required=True)
parserCore.add_argument('-submit-command', action="store", dest="submit_command", help="command used to submit job", required=True)
parserCore.add_argument('-submit-file', action="store", dest="submit_file", help="filename to a job file template, for HPC system that required a job file")

parserDev.add_argument('-nb-partitions', action="store", dest="nb_partitions", help="number of partition files per k-mer spectrums", default="0")

args =  parser.parse_args()




SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]
command = "python " + os.path.join(SCRIPT_DIR, "simka-pipeline.py")
for i in range(1, len(sys.argv)):
    if sys.argv[i][0] == "-":
        command += " " + sys.argv[i] + " "
    else:
        command += " \"" + sys.argv[i] + "\" "

#We notify that simka is run in HPC mode
command += " -hpc "

ret = os.system(command)
exit(ret)
