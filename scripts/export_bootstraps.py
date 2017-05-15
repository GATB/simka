
import os, sys, glob, shutil
#os.chdir(os.path.split(os.path.realpath(__file__))[0])


inputDir = sys.argv[1] #path to a simka "matrix_binary" dir
outputDir = sys.argv[2]
simkaExportBinFilename = sys.argv[3] #path to simka2-export binary (./scripts/simka2/bin/simka2-export)

#SIMKA2_EXPORT_BIN = "simka2/bin/simka2-export"
#command = "simka2/bin/simka2-export -out __results__/results_k31_t2/0.4/ -in __results__/results_k31_t2/matrix_binary/0.400000/"

matrixInfosFilename = os.path.join(inputDir, "matrix_infos.bin")

for dir in glob.glob(os.path.join(inputDir, "*")):
    if os.path.isdir(dir):
        #print dir
        dirname = os.path.basename(dir)
        dirname = str(float(dirname))
        #print dirname

        outputDirPath = os.path.join(outputDir, dirname)
        #print outputDirPath

        if not os.path.exists(outputDirPath):
            os.makedirs(outputDirPath, -1)

        shutil.copy(matrixInfosFilename, dir)

        command = simkaExportBinFilename + " -in " + dir + " -out " + outputDirPath
        os.system(command)
