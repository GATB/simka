
import os, argparse, shutil
from simka2_utils import Simka2ResourceAllocator



parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-kmer-size', action="store", dest="_kmerSize", default="21")
parser.add_argument('-abundance-min', action="store", dest="_abundanceMin", default="2")
parser.add_argument('-abundance-max', action="store", dest="abundance_max", help="max abundance a kmer can have to be considered", default="0")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")
parser.add_argument('-submit-command', action="store", dest="submit_command", help="command used to submit job")
parser.add_argument('-submit-file', action="store", dest="submit_file", help="filename to a job file template, for HPC system that required a job file")
parser.add_argument('-nb-partitions', action="store", dest="nb_partitions", help="number of partition files per k-mer spectrums", default="0")
parser.add_argument('-max-reads', action="store", dest="max_reads", default="0", help="maximum number of reads per sample to process")
parser.add_argument('-min-read-size', action="store", dest="min_read_size", default="0", help="minimal size a read should have to be kept")
parser.add_argument('-min-shannon-index', action="store", dest="min_read_shannon_index", default="0", help="minimal Shannon index a read should have to be kept. Float in [0,2]")
parser.add_argument('-simple-dist', action="store_true", dest="simple_dist", help="compute all simple distances (Chord, Hellinger...)")
parser.add_argument('-complex-dist', action="store_true", dest="complex_dist", help="compute all complex distances (Jensen-Shannon...)")
parser.add_argument('-nb-kmers', action="store", dest="nb_kmers", help="number of kmers used to estimate distances", default="100000")

args =  parser.parse_args()

def init_settings():

    nbPartitions = args.nb_partitions
    if nbPartitions == "0":
        resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs), None, None)
        maxJobs, coresPerJob = resourceAllocator.executeForDistanceJobs(-1)
        nbPartitions = min(200, maxJobs*coresPerJob)
    nbPartitions = max(100, nbPartitions)
    #nbPartitions = 1 #for simkaMin distance computation

    settings_file = open(os.path.join(args._databaseDir, "settings.txt"), "w")
    settings_file.write("kmer-size: " + args._kmerSize + "\n")
    settings_file.write("nb-partitions: " + str(nbPartitions)  + "\n")
    settings_file.write("abundance-min: " + args._abundanceMin  + "\n")
    if args.abundance_max == "0":
        settings_file.write("abundance-max: " + str((2**31-1))  + "\n")
    else:
        settings_file.write("abundance-max: " + args.abundance_max  + "\n")
    settings_file.write("max-reads: " + args.max_reads  + "\n")
    settings_file.write("min-read-size: " + args.min_read_size  + "\n")
    settings_file.write("min-shannon-index: " + args.min_read_shannon_index  + "\n")
    settings_file.write("simple-dist: " + ("1" if args.simple_dist else "0")  + "\n")
    settings_file.write("complex-dist: " + ("1" if args.complex_dist else "0")  + "\n")
    settings_file.write("nb-kmers: " + args.nb_kmers + "\n")

    settings_file.close()


#def init_distance_matrix():
    #distance_matrix_dir = os.path.join(args._databaseDir, "distance_matrix")
    #os.makedirs(distance_matrix_dir)

    #open(os.path.join(distance_matrix_dir, "processed_ids.bin"), "w").close()
    #os.makedirs(os.path.join(distance_matrix_dir, "mat_abundance_braycurtis"))
    #os.makedirs(os.path.join(distance_matrix_dir, "mat_presenceAbsence_jaccard"))

    #os.makedirs(os.path.join(args._databaseDir, "distance_computation"))

def main():

    settings_filename = os.path.join(args._databaseDir, "settings.txt")
    if os.path.exists(settings_filename):
        #raise Exception("Can't create database (" + args._databaseDir + "). This directory already exists.")
        exit(0)
    else:
        if(os.path.exists(args._databaseDir)):
            shutil.rmtree(args._databaseDir, ignore_errors=True)
        os.makedirs(args._databaseDir)

    init_settings()

main()