
import os, argparse, shutil
from simka2_utils import Simka2ResourceAllocator



parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-kmer-size', action="store", dest="_kmerSize", default="31")
parser.add_argument('-abundance-min', action="store", dest="_abundanceMin", default="2")
parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")
parser.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of job that can be submitted at a given time", default="0")
parser.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")
parser.add_argument('-submit-command', action="store", dest="submit_command", help="command used to submit job")
parser.add_argument('-submit-file', action="store", dest="submit_file", help="filename to a job file template, for HPC system that required a job file")

args =  parser.parse_args()


def init_settings():

    resourceAllocator = Simka2ResourceAllocator(bool(args._isHPC), int(args._nbCores), int(args._maxMemory), int(args._maxJobs), None, None)
    maxJobs, coresPerJob = resourceAllocator.executeForDistanceJobs(-1)
    nbPartitions = min(200, maxJobs)

    print maxJobs
    settings_file = open(os.path.join(args._databaseDir, "settings.txt"), "w")
    settings_file.write("kmer-size: " + args._kmerSize + "\n")
    settings_file.write("nb-partitions: " + str(nbPartitions)  + "\n")
    settings_file.write("abundance-min: " + args._abundanceMin  + "\n")
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
            shutil.rmtree(args._databaseDir)
        os.makedirs(args._databaseDir)

    init_settings()

main()