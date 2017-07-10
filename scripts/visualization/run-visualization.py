
#python create_heatmaps.py matrixFolder simkaRscriptFolder
import os
from os import listdir
from os.path import isfile, join, splitext
import sys, argparse


parser = argparse.ArgumentParser()


parserFile = parser.add_argument_group("in/out options")
parserMetadata = parser.add_argument_group("metadata options")
parserVisualization = parser.add_argument_group("visualization options")

parserFile.add_argument('-in', action="store", dest="input_dir", help="simka result directory containing distance matrices", required=True)
parserFile.add_argument('-out', action="store", dest="output_dir", help="output directory for figures", required=True)
parserFile.add_argument('-width', action="store", dest="width", help="width of figures in inches", default="7")
parserFile.add_argument('-height', action="store", dest="height", help="height of figures in inches", default="7")
parserFile.add_argument('-format', action="store", dest="format", help="output format (pdf or png)", default="png")

parserMetadata.add_argument('-metadata-in', action="store", dest="metadata_filename", help="filename containing metadata of datasets in csv format (separator=;)")
parserMetadata.add_argument('-metadata-variable', action="store", dest="metadata_variable", help="the name of the variable you want to use in metadata table (column name)")

parserVisualization.add_argument('-heatmap', action="store_true", dest="want_heatmap", help="compute and output heatmap")
parserVisualization.add_argument('-tree', action="store_true", dest="want_tree", help="compute and output hierachical clustering as a dendrogram")
parserVisualization.add_argument('-pca', action="store_true", dest="want_pca", help="compute and output pca (more precisely mds/PCoA)")
parserVisualization.add_argument('-pca-axis-1', action="store", dest="pca_axis_1", help="the number of the first axis of the PCA (only used if -pca)", default="1")
parserVisualization.add_argument('-pca-axis-2', action="store", dest="pca_axis_2", help="the number of the second axis of the PCA (only used if -pca)", default="2")

args =  parser.parse_args()



matrix = {}
#matrix = {
#"presenceAbsence_sorensen": [],
#"presenceAbsence_jaccard": [],
#"abundance_jaccard": [],
#"abundance_brayCurtis": [],
#}

def add_metadata_args(command):

	command += " " + args.width
	command += " " + args.height
	command += " " + args.format

	if args.metadata_filename != None:
		if args.metadata_variable == None:
			print("Please specify wanted variable name (column name) in metadata table with option -metadata-variable")
			exit(1)

		command += " " + args.metadata_filename + " " + args.metadata_variable

	return command

def outputHeatmap(outputFilename, matrixAsymFilename, matrixSymFilename):
	if not args.want_heatmap: return
#asymFilename = matrixAsymFilename + _outputFilenameSuffix + ".csv";
#	normFilename = matrixNormFilename + _outputFilenameSuffix + ".csv";
#	outputFilename = outputFilenamePrefix + _outputFilenameSuffix + ".png";
#print(matrixAsymFilename)
#	print(matrixNormFilename)
#	print(outputFilename)
#command = "Rscript " +  heatmap_script_filename + " " + join(mat_input_dir, matrixAsymFilename) + " " + join(mat_input_dir, matrixNormFilename) + " " + join(mat_input_dir, outputFilename)
	command = "Rscript " +  heatmap_script_filename + " " + join(args.input_dir, matrixAsymFilename) + " " + join(args.input_dir, matrixSymFilename) + " " + join(args.output_dir, outputFilename)
	command = add_metadata_args(command)

	print("\t"+command)
#print command
	os.system(command)# + " > /dev/null 2>&1  ")

def outputHclust(outputFilename, matrixNormFilename):
	if not args.want_tree: return

	command = "Rscript " +  hclust_script_filename + " " + join(args.input_dir, matrixNormFilename) + " " + join(args.output_dir, outputFilename)
	command = add_metadata_args(command)

	print("\t"+command)
	#print command
	os.system(command)# + " > /dev/null 2>&1  ")

def outputPca(outputFilename, matrixNormFilename):
	if not args.want_pca: return
	
	command = "Rscript " +  pca_script_filename + " " + join(args.input_dir, matrixNormFilename) + " " + join(args.output_dir, outputFilename) + " " + args.pca_axis_1 + " " + args.pca_axis_2
	command = add_metadata_args(command)


	print("\t"+command)
	#print(args.metadata_filename)
	#print(args.metadata_variable)
	#print command
	os.system(command)# + " > /dev/null 2>&1  ")

def execute():
	files = [ f for f in listdir(args.input_dir) if isfile(join(args.input_dir,f))]
	for filename in files:
		asym = False
		if not ".csv.gz" in filename: continue
		if "asym" in filename:
			asym = True
			asym_filename = filename
			filename = filename.replace("_asym", "")
		method_name = filename.split(".")[0]
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
		print(method_name)
		#one version of the similairty function (sym)
		if len(matrix_filenames) == 1:
			#print("lala")
			outputHeatmap("heatmap_" + method_name, matrix_filenames[0], matrix_filenames[0])
			outputHclust("hclust_" + method_name, matrix_filenames[0])
			outputPca("pca_" + method_name, matrix_filenames[0])
		#two version of the similarity function (sym and asym)
		else:
			sym = ""
			asym = ""
			for filename in matrix_filenames:
				if "asym" in filename:
					asym = filename
				else:
					sym = filename
			outputHeatmap("heatmap_" + method_name, asym, sym)
			outputHclust("hclust_" + method_name, sym)
			outputPca("pca_" + method_name, sym)



#args = sys.argv

#mat_input_dir = args[1]

#try:
#	rscript_dir = args[2]
#except:
#	rscript_dir = os.path.dirname(os.path.realpath(__file__))

rscript_dir = os.path.dirname(os.path.realpath(__file__))

heatmap_script_filename = join(rscript_dir, "heatmap.r")
hclust_script_filename = join(rscript_dir, "dendro.r")
pca_script_filename = join(rscript_dir, "pca.r")

if not args.want_heatmap and not args.want_pca and not args.want_tree:
	print("Please, choose at least one option among: -heatmap -tree -pca")
	exit(1)

if not os.path.exists(args.output_dir):
	os.makedirs(args.output_dir)

execute()



