
import os, sys, glob


allowedExtensions = ["fastq", "fasta", "fq", "fa", "fastq.gz", "fasta.gz", "fq.gz", "fa.gz"]





#-------------------------------------------------------------------
inputDir = sys.argv[2]
outputFilename = sys.argv[1]
appending = sys.argv[3]


#----------------------------------------------------------------



openType = ""
if appending == "1":
	openType = "a"
else:
	openType = "w"


simka_input_file = open(outputFilename, openType)
#inputFilenames = []
for ext in allowedExtensions:
	filenames = glob.glob(inputDir + "*." + ext)
	#print(filenames)
	#print(inputDir + "*." + ext + "*")
	#inputFilenames.append(filenames)

	for filename in filenames:
		sampleID = os.path.splitext(os.path.split(filename)[1])[0]
		
		#print(sampleID)
		#exit(1)
		simka_input_file.write(sampleID + ":" + filename + "\n")

#print(inputFilenames)

simka_input_file.close()
