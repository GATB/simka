
import os, sys, argparse, shutil


parser = argparse.ArgumentParser(description='Simka options')


parser.add_argument('-in', action="store", dest="_inputFilename", help="input file of samples. One sample per line: id1: filename1...")
parser.add_argument('-out', action="store", dest="_outputDir", default="./simka_results", help="output directory for result files (distance matrices)")
parser.add_argument('-out-tmp', action="store", dest="_outputDirTmp", help="output directory for temporary files")
parser.add_argument('-keep-tmp', action="store_true", dest="_keepTmp", help="keep temporary files. Allow to update existing run of simka")
parser.add_argument('-simka-bin', action="store", dest="_simkaBinDir", help="dir containing simka binaries")

parser.add_argument('-simple-dist', action="store", dest="_simpleDist", help="compute all simple distances (Chord, Hellinger...)")
parser.add_argument('-complex-dist', action="store", dest="_complexDist", help="compute all complex distances (Jensen-Shannon...)")

parser.add_argument('-kmer-size', action="store", dest="_kmerSize", help="size of a kmer", default="31")
parser.add_argument('-abundance-min', action="store", dest="_abundanceMin", help="min abundance a kmer need to be considered", default="2")
parser.add_argument('-abundance-max', action="store", dest="_abundanceMax", help="max abundance a kmer can have to be considered", default="999999999")
#parser.add_argument('-kmer-shannon-index', action="store", dest="_simpleDist")

parser.add_argument('-max-reads', action="store", dest="_maxReads", help="maximum number of reads per sample to process. Can be -1: use all reads. Can be 0: estimate it")
parser.add_argument('-min-read-size', action="store", dest="_minReadSize", help="minimal size a read should have to be kept")
parser.add_argument('-min-shannon-index', action="store", dest="_minReadShannonIndex", help="minimal Shannon index a read should have to be kept. Float in [0,2]")

parser.add_argument('-nb-cores', action="store", dest="_nbCores", help="number of cores", default="0")
parser.add_argument('-max-memory', action="store", dest="_maxMemory", help="max memory (MB)", default="8000")



if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

args =  parser.parse_args()

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

#-----------------------------------------------------------------------------
if not os.path.exists(args._outputDirTmp):
    os.makedirs(args._outputDirTmp)
databaseDir = os.path.join(args._outputDirTmp, "simka_database")
tmpComputationDir = os.path.join(args._outputDirTmp, "simka_tmp")
#if not os.path.exists(databaseDir):
#    os.makedirs(databaseDir)
if not os.path.exists(tmpComputationDir):
    os.makedirs(tmpComputationDir)
if not os.path.exists(args._outputDir):
    os.makedirs(args._outputDir)

# Init simka
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "simka_init.py") + \
    " -database-dir " + databaseDir + \
    " -kmer-size " + args._kmerSize + \
    " -abundance-min " + args._abundanceMin
os.system(command)

#Compute all k-mer spectrums
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "compute_kmer_spectrum_all.py") + \
    " -database-dir " + databaseDir + \
    " -out-tmp " + tmpComputationDir + \
    " -in " + args._inputFilename + \
    " -simka-bin " + args._simkaBinDir
os.system(command)

#Compute distance between k-mer spectrums
command = "python " + \
    os.path.join(SCRIPT_DIR, "core", "simka2-distance.py") + \
    " -database-dir " + databaseDir + \
    " -simka-bin " + args._simkaBinDir
os.system(command)

command = os.path.join(args._simkaBinDir, "simkaDistanceExport") + \
    " -out " + args._outputDir + \
    " -in " + os.path.join(databaseDir, "distance", "matrix_binary")
os.system(command)

#Remove tmp dir (k-mer spectrums, dist...)
if not args._keepTmp:
    shutil.rmtree(databaseDir)
    shutil.rmtree(tmpComputationDir)
