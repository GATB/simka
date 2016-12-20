
import os, shutil, argparse
from core.simka2_utils import ArgumentFormatterSimka, SimkaParser, SimkaCommand




parser = SimkaParser(formatter_class=ArgumentFormatterSimka)


parserMain = parser.add_argument_group("main options")
parserCore = parser.add_argument_group("core options")
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

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores", default="0")
parserCore.add_argument('-max-memory', action="store", dest="max_memory", help="max memory (MB)", default="8000")
parserCore.add_argument('-hpc', action="store_true", dest="_isHPC", help="compute with cluster or grid system")
parserCore.add_argument('-max-jobs', action="store", dest="_maxJobs", help="maximum number of jobs that can be submitted simultaneously", default="0")
parserCore.add_argument('-submit-command', action="store", dest="submit_command", help="command used to submit job")
parserCore.add_argument('-submit-file', action="store", dest="submit_file", help="filename to a job file template, for HPC system that required a job file")

parserDev.add_argument('-nb-partitions', action="store", dest="nb_partitions", help="number of partition files per k-mer spectrums", default="0")

args =  parser.parse_args()

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

#-----------------------------------------------------------------------------
# Create dirs that hold the simka database and temporary computation files
#-----------------------------------------------------------------------------
if not os.path.exists(args.output_dir_temp):
    os.makedirs(args.output_dir_temp)
databaseDir = os.path.join(args.output_dir_temp, "simka_database")
tmpComputationDir = os.path.join(args.output_dir_temp, "simka_tmp")
if not os.path.exists(tmpComputationDir):
    os.makedirs(tmpComputationDir)
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

#-----------------------------------------------------------------------------
# Init Simka database
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(SCRIPT_DIR, "core", "simka2-init.py")
command += " -database-dir " + databaseDir
command += " -kmer-size " + args.kmer_size
command += " -nb-cores " + args.nb_cores
command += " -max-memory " + args.max_memory
command += " -abundance-min " + args.abundance_min
command += " -abundance-max " + args.abundance_max
command += " -nb-partitions " + args.nb_partitions
command += " -max-reads " + args.max_reads
command += " -min-read-size " + args.min_read_size
command += " -min-shannon-index " + args.min_read_shannon_index
if args.simple_dist: command += " -simple-dist "
if args.complex_dist: command += " -complex-dist "
command = SimkaCommand.addHPCargs(command, args)
ret = os.system(command)
if ret != 0: exit(1)


#-----------------------------------------------------------------------------
# Compute k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(SCRIPT_DIR, "core", "simka2-count.py")
command += " -database-dir " + databaseDir
command += " -out-tmp " + tmpComputationDir
command += " -in " + args.input_filename
#command += " -simka-bin " + args.simka_bin_dir + \
command += " -nb-cores " + args.nb_cores
command += " -max-memory " + args.max_memory
command = SimkaCommand.addHPCargs(command, args)
ret = os.system(command)
if ret != 0: exit(1)

#-----------------------------------------------------------------------------
# Compute distances between k-mer spectrums
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(SCRIPT_DIR, "core", "simka2-distance.py")
command += " -database-dir " + databaseDir
command += " -nb-cores " + args.nb_cores
command += " -max-memory " + args.max_memory
#command += " -simka-bin " + args.simka_bin_dir
command = SimkaCommand.addHPCargs(command, args)
ret = os.system(command)
if ret != 0: exit(1)

#-----------------------------------------------------------------------------
# Copy distance matrix binaries in result dir
#-----------------------------------------------------------------------------
matrixBinaryFilename = os.path.join(databaseDir, "distance", "matrix_binary")
matrixBinaryFilenameDest = os.path.join(args.output_dir, "matrix_binary")
if os.path.exists(matrixBinaryFilenameDest):
    shutil.rmtree(matrixBinaryFilenameDest, ignore_errors=True)
shutil.copytree(matrixBinaryFilename, matrixBinaryFilenameDest)
matrixBinaryFilename = matrixBinaryFilenameDest

#-----------------------------------------------------------------------------
# Run the distance exporter
# The distance exporter reads matrices in binary format (-in)
# extracts the rows/columns depending on supplied datasets id (-in-ids)
# and outputs distance matrices in ascii format readable by R (-out)
#-----------------------------------------------------------------------------
command = os.path.join(SCRIPT_DIR, "bin", "simka2-export") + \
    " -out " + args.output_dir + \
    " -in " + matrixBinaryFilename
ret = os.system(command)
if ret != 0: exit(1)


#-----------------------------------------------------------------------------
# Remove tmp dir (k-mer spectrums, dist...)
# The computed k-mer spectrums and distance matrices are kept if -keep-tmp
# is set.
#-----------------------------------------------------------------------------
if not args.keep_tmp:
    shutil.rmtree(databaseDir, ignore_errors=True)
shutil.rmtree(tmpComputationDir, ignore_errors=True)


print("\n\n")
print("Results dir: " + args.output_dir)