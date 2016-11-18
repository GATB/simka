
import os, sys, argparse, shutil
from core.simka2_database import SimkaDatabase


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parserMain = parser.add_argument_group("main options")
parserCore = parser.add_argument_group("core options")
parserDistance = parser.add_argument_group("distance options")
parserKmer = parser.add_argument_group("k-mer options")
parserRead = parser.add_argument_group("read options")

parserMain.add_argument('-in', action="store", dest="input_filename", help="input file of samples. One sample per line: id1: filename1...", required=True)
parserMain.add_argument('-out-tmp', action="store", dest="output_dir_temp", help="output directory for temporary files", required=True)
parserMain.add_argument('-simka-bin', action="store", dest="simka_bin_dir", help="dir containing simka binaries", required=True)
parserMain.add_argument('-out', action="store", dest="output_dir", default="./simka_results", help="output directory for result files (distance matrices)")
parserMain.add_argument('-keep-tmp', action="store_true", dest="keep_tmp", help="keep temporary files. Allow to update existing run of simka")

parserDistance.add_argument('-simple-dist', action="store_true", dest="simple_dist", help="compute all simple distances (Chord, Hellinger...)")
parserDistance.add_argument('-complex-dist', action="store_true", dest="complex_dist", help="compute all complex distances (Jensen-Shannon...)")

parserKmer.add_argument('-kmer-size', action="store", dest="kmer_size", help="size of a kmer", default="31")
parserKmer.add_argument('-abundance-min', action="store", dest="abundance_min", help="min abundance a kmer need to be considered", default="2")
parserKmer.add_argument('-abundance-max', action="store", dest="abundance_max", help="max abundance a kmer can have to be considered", default="999999999")
#parser.add_argument('-kmer-shannon-index', action="store", dest="simple_dist")

parserRead.add_argument('-max-reads', action="store", dest="max_reads", help="maximum number of reads per sample to process. Can be -1: use all reads. Can be 0: estimate it")
parserRead.add_argument('-min-read-size', action="store", dest="min_read_size", help="minimal size a read should have to be kept")
parserRead.add_argument('-min-shannon-index', action="store", dest="min_read_shannon_index", help="minimal Shannon index a read should have to be kept. Float in [0,2]")

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores", default="0")
parserCore.add_argument('-max-memory', help="max memory (MB)", default="8000")


#args =  parser.parse_args()
#print args.max_memory
#exit(1)
#if len(sys.argv) == 1:
#    parser.print_help()
#    exit(1)

args =  parser.parse_args()

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

#-----------------------------------------------------------------------------
if not os.path.exists(args.output_dir_temp):
    os.makedirs(args.output_dir_temp)
databaseDir = os.path.join(args.output_dir_temp, "simka_database")
tmpComputationDir = os.path.join(args.output_dir_temp, "simka_tmp")
#if not os.path.exists(databaseDir):
#    os.makedirs(databaseDir)
if not os.path.exists(tmpComputationDir):
    os.makedirs(tmpComputationDir)
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)



#-----------------------------------------------------------------------------

# Init simka
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "simka2-init.py") + \
    " -database-dir " + databaseDir + \
    " -kmer-size " + args.kmer_size + \
    " -nb-cores " + args.nb_cores + \
    " -max-memory " + args.max_memory + \
    " -max-jobs " + "0" + \
    " -abundance-min " + args.abundance_min
os.system(command)


#Compute all k-mer spectrums
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "simka2-count.py") + \
    " -database-dir " + databaseDir + \
    " -out-tmp " + tmpComputationDir + \
    " -in " + args.input_filename + \
    " -simka-bin " + args.simka_bin_dir + \
    " -nb-cores " + args.nb_cores + \
    " -max-memory " + args.max_memory + \
    " -max-jobs " + "0"
os.system(command)

#Compute distance between k-mer spectrums
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "simka2-distance.py") + \
    " -database-dir " + databaseDir + \
    " -simka-bin " + args.simka_bin_dir
os.system(command)

#-----------------------------------------------------------------------------
# Write the list of available dataset ids on disk for distance exporter
# Dataset ids are retrieved from input filename (-in)
# It allows the final distance matrix to have dataset in the same order
# than supplied input filename (-in)
#-----------------------------------------------------------------------------
#datasetIdsFilename = os.path.join(tmpComputationDir, "__simka_dataset_ids.txt")
#datasetIdsFile = open(datasetIdsFilename, "w")
#inputFile = open(args.input_filename, "r")
#for line in inputFile:
#    line = line.strip()
#    if line == "": continue
#    id, filenames = line.replace(" ", "").split(":")
#    datasetIdsFile.write(id + "\n")
#datasetIdsFile.close()


#-----------------------------------------------------------------------------
# Copy distance matrix binaries in result dir
#-----------------------------------------------------------------------------
matrixBinaryFilename = os.path.join(databaseDir, "distance", "matrix_binary")
matrixBinaryFilenameDest = os.path.join(args.output_dir, "matrix_binary")
if os.path.exists(matrixBinaryFilenameDest):
    shutil.rmtree(matrixBinaryFilenameDest)
shutil.copytree(matrixBinaryFilename, matrixBinaryFilenameDest)
matrixBinaryFilename = matrixBinaryFilenameDest

#-----------------------------------------------------------------------------
# Run the distance exporter
# The distance exporter reads matrices in binary format (-in)
# extracts the rows/columns depending on supplied datasets id (-in-ids)
# and outputs distance matrices in ascii format readable by R (-out)
#-----------------------------------------------------------------------------
command = os.path.join(args.simka_bin_dir, "simka2-export") + \
    " -out " + args.output_dir + \
    " -in " + matrixBinaryFilename #+ \
    #" -in-ids " + datasetIdsFilename
os.system(command)


#Remove tmp dir (k-mer spectrums, dist...)
if not args.keep_tmp:
    shutil.rmtree(databaseDir)
shutil.rmtree(tmpComputationDir)


print("\n\n")
print("Results dir: " + args.output_dir)