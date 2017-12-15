
#python create_heatmaps.py matrixFolder simkaRscriptFolder
import os
from os import listdir
from os.path import isfile, join, splitext
import sys, argparse


#os.chdir(os.path.split(os.path.realpath(__file__))[0])




#-------------------------------------------------------------------------------------------------------------
# ArgumentFormatterSimka
#-------------------------------------------------------------------------------------------------------------
class SimkaParser(argparse.ArgumentParser):

    def error(self, message):
        print("")
        sys.stderr.write('error: %s\n' % message)
        print("")
        self.print_help()
        sys.exit(2)


class ArgumentFormatterSimka(argparse.HelpFormatter):


    #def _fill_text(self, text, width, indent):
    #    return ''.join([indent + line for line in text.splitlines(True)])
    def _split_lines(self, text, width):
        return text.splitlines()

    #remove default args layout
    def _format_args(self, action, default_metavar):
        result = ""
        return result

    #Remove "usage: ..." header
    def _format_usage(self, usage, actions, groups, prefix):
        return ""


    #Changed layout of each item
    def _get_help_string(self, action):

        text = ""

        if type(action) == argparse._StoreAction:
            text =  "(1 arg) :    " + action.help
        elif type(action) == argparse._StoreTrueAction:
            text =  "(0 arg) :    " + action.help

        if type(action) == argparse._StoreAction and action.default != None:
            text += " [Default: " + str(action.default) + "]"
        #print type(action), action
        #print action
        #return "-5-"
        #return action.help
        if text != "":
            return text

        return "__none__"

    #Hack for removing useless "optional arguments:" section
    def _join_parts(self, part_strings):
        #print part_strings
        return ''.join([part
                        for part in part_strings
                        if part and part is not argparse.SUPPRESS and not "optional arguments:" in part and not "__none__" in part and not "--help" in part])






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
sketchDir = os.path.join(args.out, "sketch")
if not os.path.exists(sketchDir): os.makedirs(sketchDir)
sketchFilename = os.path.join(sketchDir, "sketch.bin")
distanceOutputDir = os.path.join(args.out, "distance")
if not os.path.exists(distanceOutputDir): os.makedirs(distanceOutputDir)

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

distanceCommand = args.bin + " distance "
distanceCommand += " -in1 " + sketchFilename
distanceCommand += " -in2 " + sketchFilename
distanceCommand += " -out " + distanceOutputDir
distanceCommand += " -nb-cores " + args.nb_cores

exportCommand = args.bin + " export "
exportCommand += " -in " + distanceOutputDir
exportCommand += " -in1 " + sketchFilename
exportCommand += " -in2 " + sketchFilename
#exportCommand += " -in-ids " + distanceOutputDir #not applicable here
exportCommand += " -out " + args.out
exportCommand += " -nb-cores " + args.nb_cores


print("\n\n#-----------------------------")
print("# Sketching")
print("#-----------------------------\n")
ret = os.system(sketchCommand)
if ret != 0: print("ERROR"); exit(1)

print("\n\n#-----------------------------")
print("# Computing distances")
print("#-----------------------------\n")
ret = os.system(distanceCommand)
if ret != 0: print("ERROR"); exit(1)


print("\n\n#-----------------------------")
print("# Exporting distances")
print("#-----------------------------\n")
ret = os.system(exportCommand)
if ret != 0: print("ERROR"); exit(1)

print("\n\n")
print("Result dir: " + args.out)