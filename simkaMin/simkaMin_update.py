#!/usr/bin/env python

#*****************************************************************************
#   SimkaMin: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2019  INRIA
#   Authors: G.Benoit, C.Lemaitre, P.Peterlongo
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

import os, struct, shutil
from os import listdir
from os.path import isfile, join, splitext
import sys, argparse

from simkaMin_utils import SimkaParser, ArgumentFormatterSimka, read_sketch_header, is_executable



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

parserMain.add_argument('-in', action="store", dest="input_filename", help="input file of datasets (datasets to add to existing simka results", required=True)
parserMain.add_argument('-in-to-update', action="store", dest="input_existingResults", help="path to existing simka results to update (existing results will be overwritten)", required=True)
parserMain.add_argument('-bin', action="store", dest="bin", help="path to simkaMinCore program (to be specified if not in PATH, or not in standard installation directory <simkaDir>/build/bin/simkaMinCore)")

parserKmer.add_argument('-filter', action="store_true", dest="filter", help="filter out k-mer seen one time (potentially erroneous)")


parserRead.add_argument('-max-reads', action="store", dest="max_reads", default="0", help="maximum number of reads per sample to process")
parserRead.add_argument('-min-read-size', action="store", dest="min_read_size", default="0", help="minimal size a read should have to be kept")
parserRead.add_argument('-min-shannon-index', action="store", dest="min_shannon_index", default="0", help="minimal Shannon index a read should have to be kept. Float in [0,2]")

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores", default="0")
parserCore.add_argument('-max-memory', action="store", dest="max_memory", help="max memory (MB)", default="8000")


args =  parser.parse_args()


# Check SimkaMinCore executable
# -----------------------------

simkaMinCoreBin=args.bin
if args.bin is not None:
    # given by the user
    if not is_executable(simkaMinCoreBin):
        print("Error: "+simkaMinCoreBin+" not found or not executable, should be <SimkaDirectory>/build/bin/simkaMinCore")
        exit(1)
else:
    # Check if is in the PATH
    simkaMinCoreBin="simkaMinCore"
    if not is_executable(simkaMinCoreBin):
        # not in PATH, checking "../build/bin/simkaMinCore"
        simkaMinCoreBin=os.path.join(os.path.split(os.path.realpath(__file__))[0],"../build/bin/simkaMinCore")
        if not is_executable(simkaMinCoreBin):
            print("Error: simkaMinCore executable not found, please give the executable path with option -bin (should be <SimkaDirectory>/build/bin/simkaMinCore)")
            exit(1)



#-------------------------------------------------------------------------------------------------------------
# SimkaMin pipeline
#-------------------------------------------------------------------------------------------------------------

#Create some dirs and filenames
#if not os.path.exists(args.out): os.makedirs(args.out)
existingDir = args.input_existingResults
sketchDir = os.path.join(existingDir, "sketch")
#if not os.path.exists(sketchDir): os.makedirs(sketchDir)
sketchFilename_existing = os.path.join(sketchDir, "sketch.bin")
sketchFilename_new = os.path.join(sketchDir, "sketch_new.bin")
distanceOutputDir = os.path.join(existingDir, "distance")

distanceDir_existingVsNew = os.path.join(distanceOutputDir, "existingVsNew")
if not os.path.exists(distanceDir_existingVsNew): os.makedirs(distanceDir_existingVsNew)
distanceDir_newVsNew = os.path.join(distanceOutputDir, "newVsNew")
if not os.path.exists(distanceDir_newVsNew): os.makedirs(distanceDir_newVsNew)


#simkaMin_pipeline_filename = "./simkaMin_pipeline.py"



#Existing datasets: datasets that have already been processed by SimkaMin
#   - ids and k-mers are contained in (-in-existing)/sketch/sketch.bin
#   - distances are contained in (-in-existing)/distance/mat_*.bin
#New datasets: datasets to add contains in -in file


existing_sketch_header = read_sketch_header(sketchFilename_existing)
print(existing_sketch_header)



#Sketch new datasets
command_sketchNewDatasets = simkaMinCoreBin + " sketch "
command_sketchNewDatasets += " -in " + args.input_filename
command_sketchNewDatasets += " -out " + sketchFilename_new
command_sketchNewDatasets += " -seed " + str(existing_sketch_header["seed"])
command_sketchNewDatasets += " -kmer-size " + str(existing_sketch_header["kmerSize"])
command_sketchNewDatasets += " -nb-kmers " + str(existing_sketch_header["sketchSize"])
if args.filter: command_sketchNewDatasets += " -filter "
command_sketchNewDatasets += " -max-reads " + args.max_reads
command_sketchNewDatasets += " -min-read-size " + args.min_read_size
command_sketchNewDatasets += " -min-shannon-index " + args.min_shannon_index
command_sketchNewDatasets += " -nb-cores " + args.nb_cores
command_sketchNewDatasets += " -max-memory " + args.max_memory

#Compute distance between existing datasets and new datasets
command_distance_existingVsNew = simkaMinCoreBin + " distance "
command_distance_existingVsNew += " -in1 " + sketchFilename_existing
command_distance_existingVsNew += " -in2 " + sketchFilename_new
command_distance_existingVsNew += " -out " + distanceDir_existingVsNew
command_distance_existingVsNew += " -nb-cores " + args.nb_cores



#Compute distance between new datasets and new datasets
command_distance_newVsNew = simkaMinCoreBin + " distance "
command_distance_newVsNew += " -in1 " + sketchFilename_new
command_distance_newVsNew += " -in2 " + sketchFilename_new
command_distance_newVsNew += " -out " + distanceDir_newVsNew
command_distance_newVsNew += " -nb-cores " + args.nb_cores

#Update existing distance matrix
command_distanceMatrix_update = simkaMinCoreBin + " matrix-update "
command_distanceMatrix_update += " -in " + distanceOutputDir
command_distanceMatrix_update += " -in1 " + sketchFilename_existing
command_distanceMatrix_update += " -in2 " + sketchFilename_new

#Append new sketch to existing sketch
command_sketch_append = simkaMinCoreBin + " append "
command_sketch_append += " -in1 " + sketchFilename_existing
command_sketch_append += " -in2 " + sketchFilename_new


exportCommand = simkaMinCoreBin + " export "
exportCommand += " -in " + distanceOutputDir
exportCommand += " -in1 " + sketchFilename_existing
exportCommand += " -in2 " + sketchFilename_existing
#exportCommand += " -in-ids " + distanceOutputDir #not applicable here
exportCommand += " -out " + args.input_existingResults


print("\n\n#-----------------------------")
print("# Sketching new datasets")
print("#-----------------------------\n")
ret = os.system(command_sketchNewDatasets)
if ret != 0: print("ERROR"); exit(1)

print("\n\n#-----------------------------")
print("# Computing distances between existing datasets and new datasets")
print("#-----------------------------\n")
ret = os.system(command_distance_existingVsNew)
if ret != 0: print("ERROR"); exit(1)

########################
#exportCommand = args.bin + " export "
#exportCommand += " -in " + distanceDir_existingVsNew
#exportCommand += " -in1 " + sketchFilename_existing
#exportCommand += " -in2 " + sketchFilename_new
#exportCommand += " -in-ids " + distanceOutputDir #not applicable here
#exportCommand += " -out " + distanceDir_existingVsNew

#os.system(exportCommand)
#os.system("gzip -cd "+ distanceDir_existingVsNew +"/mat_abundance_braycurtis.csv.gz")
########################

print("\n\n#-----------------------------")
print("# Computing distances between new datasets")
print("#-----------------------------\n")
ret = os.system(command_distance_newVsNew)
if ret != 0: print("ERROR"); exit(1)

print("\n\n#-----------------------------")
print("# Update existing distance matrices")
print("#-----------------------------\n")
ret = os.system(command_distanceMatrix_update)
if ret != 0: print("ERROR"); exit(1)

print("\n\n#-----------------------------")
print("# Append new sketch to existing sketch")
print("#-----------------------------\n")
ret = os.system(command_sketch_append)
if ret != 0: print("ERROR"); exit(1)

print("\n\n#-----------------------------")
print("# Exporting distances")
print("#-----------------------------\n")
ret = os.system(exportCommand)
if ret != 0: print("ERROR"); exit(1)

#Clear temp dir
shutil.rmtree(distanceDir_existingVsNew)
shutil.rmtree(distanceDir_newVsNew)
os.remove(sketchFilename_new)

print("\n\n")
print("Result dir: " + existingDir)
