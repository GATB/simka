
import os, argparse

parser = argparse.ArgumentParser(description='Description')

parser.add_argument('-database-dir', action="store", dest="_databaseDir")
parser.add_argument('-kmer-size', action="store", dest="_kmerSize", default="31")
parser.add_argument('-abundance-min', action="store", dest="_abundanceMin", default="2")

args =  parser.parse_args()


def init_settings():
    settings_file = open(os.path.join(args._databaseDir, "settings.txt"), "w")
    settings_file.write("kmer-size: " + args._kmerSize + "\n")
    settings_file.write("nb-partitions: " + "200"  + "\n")
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
    if os.path.exists(args._databaseDir):
        #raise Exception("Can't create database (" + args._databaseDir + "). This directory already exists.")
        exit(0)
    else:
        os.makedirs(args._databaseDir)

    init_settings()
    #init_distance_matrix()


main()