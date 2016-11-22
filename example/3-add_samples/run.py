
import os, sys, random, shutil

scriptDir = sys.argv[1]
binDir = sys.argv[2]
inputFilename = sys.argv[3]
maxSamples = max(int(sys.argv[4]), 1)


CURRENT_DIR = os.path.split(os.path.realpath(__file__))[0]

#----------------------------------------------
def clear():
	filename = os.path.join(os.path.split(inputFilename)[0], "______next_simka_input.txt")
	if os.path.exists(filename): os.remove(filename)
	
	#dir = os.path.join(CURRENT_DIR, "simka_results")
	#if os.path.exists(dir): shutil.rmtree(dir)
	dir = os.path.join(CURRENT_DIR, "simka_tmp")
	if os.path.exists(dir): shutil.rmtree(dir)

clear()
#----------------------------------------------

nbSamplesToProcess = 0

inputFile = open(inputFilename, "r")
for line in inputFile:
	line = line.strip()
	if line == "": continue

	nbSamplesToProcess += 1



#----------------------------------------------
def createNextInput(nbSamplesToPick, passId):

	#filename = os.path.join(CURRENT_DIR, str(passId) + "__next_simka_input.txt")
	filename = os.path.join(os.path.split(inputFilename)[0], "______next_simka_input.txt")
	if os.path.exists(filename):
		os.remove(filename)
	
	file = open(filename, "w")
	
	nbSamplesPicked = 0
	
	for line in inputFile:
		line = line.strip()
		if line == "": continue
	
		file.write(line + "\n")
		nbSamplesPicked += 1

		if nbSamplesPicked >= nbSamplesToPick: break

	file.close()

	return filename


inputFile.seek(0)
nbSamplesProcessed = 0
passId = 0

while(nbSamplesProcessed < nbSamplesToProcess):
	nbSamples = random.randint(1, maxSamples)
	
	nbSamplesRemaining = nbSamplesToProcess - nbSamplesProcessed
	nbSamples = min(nbSamples, nbSamplesRemaining)
	
	filename = createNextInput(nbSamples, passId)
	
	command = "python " + os.path.join(scriptDir, "simka.py") + \
		" -in " + filename + \
		" -out-tmp " + os.path.join(CURRENT_DIR, "simka_tmp") + \
		" -simka-bin " + binDir + \
		" -out " + os.path.join(CURRENT_DIR, "simka_results") + \
		" -keep-tmp "

	os.system(command)
	#print command

	#print nbSamples

	nbSamplesProcessed += nbSamples

	passId += 1



inputFile.close()
clear()