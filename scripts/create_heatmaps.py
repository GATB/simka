
#python create_heatmaps.py matrixFolder simkaRscriptFolder
import os
from os import listdir
from os.path import isfile, join, splitext
import sys

matrix = {}
#matrix = {
#"presenceAbsence_sorensen": [],
#"presenceAbsence_jaccard": [],
#"abundance_jaccard": [],
#"abundance_brayCurtis": [],
#}


def outputHeatmap(outputFilename, matrixAsymFilename, matrixSymFilename):

#asymFilename = matrixAsymFilename + _outputFilenameSuffix + ".csv";
#	normFilename = matrixNormFilename + _outputFilenameSuffix + ".csv";
#	outputFilename = outputFilenamePrefix + _outputFilenameSuffix + ".png";
#print(matrixAsymFilename)
#	print(matrixNormFilename)
#	print(outputFilename)
#command = "Rscript " +  heatmap_script_filename + " " + join(mat_input_dir, matrixAsymFilename) + " " + join(mat_input_dir, matrixNormFilename) + " " + join(mat_input_dir, outputFilename)
	command = "Rscript " +  heatmap_script_filename + " " + join(mat_input_dir, matrixAsymFilename) + " " + join(mat_input_dir, matrixSymFilename) + " " + join(mat_input_dir, outputFilename)
	print(command)
#print command
	os.system(command)

def outputHclust(outputFilename, matrixNormFilename):
	
	command = "Rscript " +  hclust_script_filename + " " + join(mat_input_dir, matrixNormFilename) + " " + join(mat_input_dir, outputFilename)
	print(command)
	#print command
	os.system(command)

def execute():
	files = [ f for f in listdir(mat_input_dir) if isfile(join(mat_input_dir,f))]
	for filename in files:
		asym = False
		if not ".csv" in filename: continue
		if "asym" in filename:
			asym = True
			asym_filename = filename
			filename = filename.replace("_asym", "")
		method_name = splitext(filename)[0]
		method_name = method_name.replace("mat_", "")
		try:
			if asym:
				matrix[method_name].append(asym_filename)
			else:
				matrix[method_name].append(filename)
		except:
			matrix[method_name] = []
			if asym:
				matrix[method_name].append(asym_filename)
			else:
				matrix[method_name].append(filename)
		#for method_name in matrix.keys():
			#print(filename, method_name)
			#if method_name in filename:
				#matrix[method_name].append(filename)
				#break

	for method_name, matrix_filenames in matrix.items():
		print("")
		#one version of the similairty function (sym)
		if len(matrix_filenames) == 1:
			#print("lala")
			outputHeatmap("heatmap_" + method_name + ".pdf", matrix_filenames[0], matrix_filenames[0])
			outputHclust("hclust_" + method_name + ".pdf", matrix_filenames[0])
		#two version of the similarity function (sym and asym)
		else:
			sym = ""
			asym = ""
			for filename in matrix_filenames:
				if "asym" in filename:
					asym = filename
				else:
					sym = filename
			outputHeatmap("heatmap_" + method_name + ".pdf", asym, sym)
			outputHclust("hclust_" + method_name + ".pdf", sym)


#outputHeatmap()
	#__outputHeatmap("heatmap_presenceAbsence_sorensen", "mat_presenceAbsence_sorensen", "mat_presenceAbsence_sorensen");
	#__outputHeatmap("heatmap_presenceAbsence_jaccard", "mat_presenceAbsence_jaccard_asym", "mat_presenceAbsence_jaccard");
	#__outputHeatmap("heatmap_abundance_jaccard", "mat_abundance_jaccard_asym", "mat_abundance_jaccard");
	#__outputHeatmap("heatmap_abundance_brayCurtis", "mat_abundance_brayCurtis", "mat_abundance_brayCurtis");


args = sys.argv

mat_input_dir = args[1]


try:
	rscript_dir = args[2]
except:
	rscript_dir = os.path.dirname(os.path.realpath(__file__))

heatmap_script_filename = join(rscript_dir, "heatmap.r")
hclust_script_filename = join(rscript_dir, "dendro.r")

execute()