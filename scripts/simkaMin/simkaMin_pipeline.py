
#python create_heatmaps.py matrixFolder simkaRscriptFolder
import os, math
from os import listdir
from os.path import isfile, join, splitext
import sys, argparse

from simkaMin_utils import SimkaParser, ArgumentFormatterSimka, read_sketch_header, ProgressBar

#os.chdir(os.path.split(os.path.realpath(__file__))[0])









#-------------------------------------------------------------------------------------------------------------
# Arg parser
#-------------------------------------------------------------------------------------------------------------
parser = SimkaParser(formatter_class=ArgumentFormatterSimka)


parserMain = parser.add_argument_group("[main options]")
parserCore = parser.add_argument_group("[core options]")
parserDistance = parser.add_argument_group("[distance options]")
parserKmer = parser.add_argument_group("[k-mer options]")
parserRead = parser.add_argument_group("[read options]")
parserDev = parser.add_argument_group("[advanced (developer) options]")

parserMain.add_argument('-in', action="store", dest="input_filename", help="input file of datasets. One sample per line: id1: filename1...", required=True)
parserMain.add_argument('-bin', action="store", dest="bin", help="path to simkaMin program (should be at build/bin/simkaMin)", required=True)
parserMain.add_argument('-out', action="store", dest="out", default="./simka_results", help="output directory for result files (distance matrices)")
parserMain.add_argument('-seed', action="store", dest="seed", default="100", help="seed used for random k-mer selection")

parserKmer.add_argument('-kmer-size', action="store", dest="kmer_size", help="size of a kmer", default="21")
parserKmer.add_argument('-nb-kmers', action="store", dest="nb_kmers", help="number of kmers used to compute distances", default="100000")
parserKmer.add_argument('-filter', action="store_true", dest="filter", help="filter out k-mer seen one time (potentially erroneous)")


parserRead.add_argument('-max-reads', action="store", dest="max_reads", default="0", help="maximum number of reads per sample to process")
parserRead.add_argument('-min-read-size', action="store", dest="min_read_size", default="0", help="minimal size a read should have to be kept")
parserRead.add_argument('-min-shannon-index', action="store", dest="min_shannon_index", default="0", help="minimal Shannon index a read should have to be kept. Float in [0,2]")

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores", default="0")
parserCore.add_argument('-max-memory', action="store", dest="max_memory", help="max memory (MB)", default="8000")


args =  parser.parse_args()





#-------------------------------------------------------------------------------------------------------------
# SimkaMin pipeline
#-------------------------------------------------------------------------------------------------------------

#Create some dirs and filenames
if not os.path.exists(args.out): os.makedirs(args.out)
outDir = os.path.join(args.out, "simkamin")
if not os.path.exists(outDir): os.makedirs(outDir)
sketchDir = os.path.join(outDir, "sketch")
if not os.path.exists(sketchDir): os.makedirs(sketchDir)
sketchFilename = os.path.join(sketchDir, "sketch.bin")
distanceOutputDir = os.path.join(outDir, "distance")
if not os.path.exists(distanceOutputDir): os.makedirs(distanceOutputDir)
logsDir = os.path.join(outDir, "logs")
if not os.path.exists(logsDir): os.makedirs(logsDir)

#Create commands
sketchCommand = args.bin + " sketch "
sketchCommand += " -in " + args.input_filename
sketchCommand += " -out " + sketchFilename
sketchCommand += " -seed " + args.seed
sketchCommand += " -kmer-size " + args.kmer_size
sketchCommand += " -nb-kmers " + args.nb_kmers
if args.filter: sketchCommand += " -filter "
sketchCommand += " -max-reads " + args.max_reads
sketchCommand += " -min-read-size " + args.min_read_size
sketchCommand += " -min-shannon-index " + args.min_shannon_index
sketchCommand += " -nb-cores " + args.nb_cores
sketchCommand += " -max-memory " + args.max_memory



exportCommand = args.bin + " export "
exportCommand += " -in " + distanceOutputDir
exportCommand += " -in1 " + sketchFilename
exportCommand += " -in2 " + sketchFilename
#exportCommand += " -in-ids " + distanceOutputDir #not applicable here
exportCommand += " -out " + args.out
exportCommand += " -nb-cores " + args.nb_cores


#print("\n\n#-----------------------------")
#print("# Sketching")
#print("#-----------------------------\n")
print("\n\n")
ret = os.system(sketchCommand)
if ret != 0: print("ERROR"); exit(1)











#print("\n\n#-----------------------------")
#print("# Computing distances")
#print("#-----------------------------\n")
print("\n\n")

#Create binary matrix file (required in case the following distance commands are run in parallel
if os.path.exists(distanceOutputDir + "/mat_presenceAbsence_jaccard.bin"): os.remove(distanceOutputDir + "/mat_presenceAbsence_jaccard.bin")
if os.path.exists(distanceOutputDir + "/mat_abundance_braycurtis.bin"): os.remove(distanceOutputDir + "/mat_abundance_braycurtis.bin")
open(distanceOutputDir + "/mat_presenceAbsence_jaccard.bin", "wb").close()
open(distanceOutputDir + "/mat_abundance_braycurtis.bin", "wb").close()

sketch_header = read_sketch_header(sketchFilename)
nbDatasetToProcess = sketch_header["nbDatasets"]
MAX_DATASETS_PROCESS = 100
def create_distance_command(i, j, n1, n2):
    distanceCommand = args.bin + " distance "
    distanceCommand += " -in1 " + sketchFilename
    distanceCommand += " -in2 " + sketchFilename
    distanceCommand += " -out " + distanceOutputDir
    distanceCommand += " -nb-cores " + args.nb_cores
    distanceCommand += " -start-i " + str(i*MAX_DATASETS_PROCESS)
    distanceCommand += " -start-j " + str(j*MAX_DATASETS_PROCESS)
    distanceCommand += " -n-i " + str(n1)
    distanceCommand += " -n-j " + str(n2)
    distanceCommand += " > " + os.path.join(logsDir, "log_distance_" + str(i) + "-" + str(j)) + " 2>&1 "
    return distanceCommand



step = int(math.ceil( float(nbDatasetToProcess) / float(MAX_DATASETS_PROCESS)))
nbCommands = int(math.ceil( float(step * step) / float(2)))
progressBar = ProgressBar("Computing distances", nbCommands)
progressBar.start()
done = False
for i in range(0, step):
    n1 = min(MAX_DATASETS_PROCESS, nbDatasetToProcess-i*MAX_DATASETS_PROCESS)
    for j in range(i, step):
        n2 = min(MAX_DATASETS_PROCESS, nbDatasetToProcess-j*MAX_DATASETS_PROCESS)
        distanceCommand = create_distance_command(i, j, n1, n2)
        #print distanceCommand
        ret = os.system(distanceCommand)
        if ret != 0: print("ERROR"); exit(1)
        progressBar.step(1)





#print("\n\n#-----------------------------")
#print("# Exporting distances")
#print("#-----------------------------\n")
print("\n\nExporting distance matrices in csv.gz format...")
ret = os.system(exportCommand)
if ret != 0: print("ERROR"); exit(1)

print("\n\n")
print("Result dir: " + args.out)