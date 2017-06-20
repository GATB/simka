

import os, sys
from os.path import join, isfile, splitext
import shutil
from os import listdir
import argparse
import shlex



args = sys.argv

nbBoostrap = int(args[1])
outputDir = args[2]
simkaCommand = args[3]

splitCommand = shlex.split(simkaCommand)
print splitCommand
simkaOutputDir = ""
outFieldId = -1
for i in range(0, len(splitCommand)):
    field = splitCommand[i]
    if field == "-out":
        outFieldId = i+1
        simkaOutputDir = splitCommand[outFieldId]
        if not os.path.exists(simkaOutputDir):
            os.mkdir(simkaOutputDir, -1)
        simkaOutputDir = join(simkaOutputDir, "simkaBootstrap")
        if os.path.exists(simkaOutputDir):
            shutil.rmtree(simkaOutputDir)
        os.mkdir(simkaOutputDir, -1)
        splitCommand[outFieldId] = simkaOutputDir
        break
if simkaOutputDir == "":
    print("-out argument of simka is missing")
    sys.exit(0)
simkaCommand = ""
for field in splitCommand:
    simkaCommand += field + " "

if not os.path.exists(outputDir):
    os.mkdir(outputDir, -1)

def create_dir(dir_name):
    path = join(outputDir, dir_name)
    if not os.path.exists(path):
        os.mkdir(path, -1)
    #shutil.rmtree(path)
    return path

def getBoostrapOffset(distanceDir):
    ids = []
    filenames = listdir(distanceDir)

    if len(filenames) == 0:
        return 0

    for filename in filenames:
        basename, ext = splitext(filename)
        fields = basename.split("_")
        id = int(fields[-1])
        ids.append(id)
    ids.sort()
    return ids[-1] + 1


idOffset = -1
for i in range(0, nbBoostrap):

    command = simkaCommand
    os.system(command)

    for filename in listdir(simkaOutputDir):
        completeFilename = join(simkaOutputDir, filename)
        #print(filename)
        distanceName = filename.replace("mat_", "")
        distanceName = distanceName.replace(".csv", "")
        #print(distanceName)

        distanceDir = create_dir(distanceName)

        if idOffset == -1:
            idOffset = getBoostrapOffset(distanceDir)
        id = idOffset + i

        destFilename = join(distanceDir, "mat_" + distanceName + "_" + str(id) + ".csv")

        shutil.move(completeFilename, destFilename)



