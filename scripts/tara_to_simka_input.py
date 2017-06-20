"""
python tara_to_simka_input.py filenames_file metadata_file output_filename
"""

from itertools import izip
import sys
import re

args = sys.argv

filenames = args[1]
metadata = args[2]
outputFilename = args[3]

metadataFile = open(metadata).readlines()
filenameFile = open(filenames).readlines()
outputFile = open(outputFilename, "w")

filenameLines = []
for filenameLine in filenameFile:
	filenameLines.append(filenameLine.strip())

foundFilenames = {}
output_lines = []


def natural_sort(l):
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
	return sorted(l, key = alphanum_key)

def getId(id):
	if "DCM" in id:
		l = id.split("DCM")
		return l[0] + "_DCM"
	elif "SUR" in id:
		l = id.split("SUR")
		return l[0] + "_SUR"
	raise Exception("getId(" + id + "), DCM or SUR not found in id")

#/ccc/scratch/cont007/fg0001/vanniet/compareads_All_0.8-5/fastq/AHX_AAAOSU_7_1_C2EV3ACXX_clean.fastq---7DCM1GGMM11---201478680
for line in metadataFile:
	line = line.strip()
	lineSplit = line.split("---")
	id = lineSplit[1]
	#print(lineSplit)
	filenameSplit = lineSplit[0].split("/")
	#print(filenameSplit)
	filename = filenameSplit[len(filenameSplit)-1]
	
	#print(lineSplit)
	#print(filename)

	for filenameLine in filenameLines:
		if filename in filenameLine:
			foundFilenames[filenameLine] = True
			s = getId(id) + " " + filenameLine + "\n"
			#print(s)
			output_lines.append(s)
			break

output_lines_sorted = natural_sort(output_lines)
for line in output_lines_sorted:
	outputFile.write(line)

print("File not found:")
for filenameLine in filenameLines:
	if not (filenameLine in foundFilenames):
		print("\t" + filenameLine)


outputFile.close()
"""
with open(filenames) as filenamesFile
	for x, y in izip(filenamesFile, metadataFile):
		x = x.strip()
		y = y.strip()
		print("{0}\t{1}".format(x, y))
"""