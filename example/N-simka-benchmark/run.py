
import os, sys, random, shutil, time

scriptDir = sys.argv[1]
binDir = sys.argv[2]
inputFilename = sys.argv[3]
#nbSamplePerPass = max(int(sys.argv[4]), 1)
#nbSamplePerPassVariance = int(sys.argv[5])

NB_SAMPLES_PER_PASS = [1, 2, 5, 10, 50, 100, 200, 400, 800, 1600, 3200]







CURRENT_DIR = os.path.split(os.path.realpath(__file__))[0]

#----------------------------------------------
def clear():
	filename = os.path.join(os.path.split(inputFilename)[0], "______next_simka_input.txt")
	if os.path.exists(filename): os.remove(filename)

	dir = os.path.join(CURRENT_DIR, "simka_results")
	if os.path.exists(dir): shutil.rmtree(dir)

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

def write_bench(results):
	benchTotalTimeFilename = os.path.join(CURRENT_DIR, "bench_total_time.csv")
	benchPassTimeFilename = os.path.join(CURRENT_DIR, "bench_pass_time.csv")
	benchTotalTimeFile = open(benchTotalTimeFilename, "w")
	benchPassTimeFile = open(benchPassTimeFilename, "w")
	benchTotalTimeFile.write(header + "\n")
	benchPassTimeFile.write(header + "\n")
	
	for nbSamplesProcessed in sorted(results.keys()):
		r = results[nbSamplesProcessed]
		line1 = str(nbSamplesProcessed)
		line2 = str(nbSamplesProcessed)
		for nbSamplePerPass in NB_SAMPLES_PER_PASS:
			if nbSamplePerPass in r:
				line1 += ";" + str(r[nbSamplePerPass][0])
				line2 += ";" + str(r[nbSamplePerPass][1])
			else:
				line1 += ";"
				line2 += ";"
		benchPassTimeFile.write(line1+"\n")
		benchTotalTimeFile.write(line2+"\n")
	
	benchTotalTimeFile.close()
	benchPassTimeFile.close()


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

results = {}

header = "nb_datasets"
for nbSamplePerPass in NB_SAMPLES_PER_PASS:

	header += ";" + "add_" + str(nbSamplePerPass)

	#if nbSamplePerPass > nbSamplesToProcess: continue
	clear()
	inputFile.seek(0)
	nbSamplesProcessed = 0
	passId = 0

	totalTimeStart = time.time()

	while(nbSamplesProcessed < nbSamplesToProcess):
		nbSamples = nbSamplePerPass

		nbSamplesRemaining = nbSamplesToProcess - nbSamplesProcessed
		if nbSamplesRemaining < nbSamplePerPass: break

		nbSamples = min(nbSamples, nbSamplesRemaining)

		filename = createNextInput(nbSamples, passId)
		passTimeStart = time.time()
		command = "python " + os.path.join(scriptDir, "simka.py") + \
			" -in " + filename + \
			" -out-tmp " + os.path.join(CURRENT_DIR, "simka_tmp") + \
			" -simka-bin " + binDir + \
			" -out " + os.path.join(CURRENT_DIR, "simka_results") + \
			" -nb-cores 1 " + \
			" -keep-tmp > /dev/null 2>&1 "

		os.system(command)
		passTime = time.time() - passTimeStart
		totalTime = time.time() - totalTimeStart
		#print command

		#print nbSamples

		nbSamplesProcessed += nbSamples

		#benchFile.write(str(nbSamplesProcessed) + ";" + str(passTime) + ";" + str(totalTime) + "\n")
		#benchFile.flush()

		if not nbSamplesProcessed in results:
			results[nbSamplesProcessed] = {}
		results[nbSamplesProcessed][nbSamplePerPass] = (passTime/float(nbSamplePerPass), totalTime)

		passId += 1

	write_bench(results)



inputFile.close()
clear()