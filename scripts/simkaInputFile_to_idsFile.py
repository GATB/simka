
#usage:
#python path/to/simka_input.txt path/to/ids.txt





import os, sys
from os.path import join, isfile

inputFilename = sys.argv[1]
outputFilename = sys.argv[2]


datasetIdsFile = open(outputFilename, "w")
inputFile = open(inputFilename, "r")

for line in inputFile:
    line = line.strip()
    if line == "": continue
    id, filenames = line.replace(" ", "").split(":")
    datasetIdsFile.write(id + "\n")

datasetIdsFile.close()