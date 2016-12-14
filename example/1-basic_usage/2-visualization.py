

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])


print("\nA collection of scripts is available to visualize simka results with R")
print("\tThose scripts are located in \"scripts/visualization\" directory")
print("\tThey can all be automatically run on simka results using the script run-visualization.py")

#-----------------------------------------------------------------------------
# Set required simka parameters
#-----------------------------------------------------------------------------
visualization_dir = "../../scripts/visualization"
distance_matrices_dir = "simka_results"
output_dir_without_color = "visualization_without_color"
output_dir_with_color = "visualization_with_color"

if not os.path.exists(distance_matrices_dir):
    print("Please run first example before: 1-simple_test.py")
    exit(1)

#-----------------------------------------------------------------------------
# Create the command to execute
#-----------------------------------------------------------------------------
command = "python "
command += visualization_dir + "/run-visualization.py"
command += " -in " + distance_matrices_dir
command += " -out " + output_dir_without_color
command += " -pca "
command += " -heatmap "
command += " -tree "
os.system(command + " > /dev/null 2>&1  ")

print("\n\n--------------------------------------------------------------")
print("- 1) Basic visualization")
print("--------------------------------------------------------------")
print("\nYou can now visualize simka results as heatmap, hierarchical clustering and pca in directory: " + output_dir_without_color)
print("\tcommand used:")
print("\t" + command)

print("\nYou have to choose at least one option among -pca  -heatmap  -tree")
print("\t-pca: compute and output PCA (more precisely MDS/PCoA)")
print("\t-heatmap: compute and output heatmap")
print("\t-tree: compute and output hierachical clustering as a dendrogram")

print("\nNote: If you use -pca, you can change pca axis using options -pca-axis-1 and -pca-axis-2")



#-----------------------------------------------------------------------------
# Run visualization script with color annotations
#-----------------------------------------------------------------------------
command = "python "
command += visualization_dir + "/run-visualization.py"
command += " -in " + distance_matrices_dir
command += " -out " + output_dir_with_color
command += " -pca "
command += " -heatmap "
command += " -tree "
command += " -metadata-in ../data/dataset_metadata.csv"
command += " -metadata-variable VARIABLE_1"
os.system(command + " > /dev/null 2>&1  ")


print("\n\n--------------------------------------------------------------")
print("- 2) Advanced visualization using metadata")
print("--------------------------------------------------------------")
print("\nYou may want to add annotations to your figures (i.e. color corresponding to given metadata)")
print("Visualization scripts accept metadata table in the given format:")
print("\n\tDATASET_ID;VARIABLE_NAME_1;VARIABLE_NAME_2")
print("\tA;1;aquatic")
print("\tB;1;human")
print("\tC;2;human")
print("\tD;2;soil")
print("\tE;3;soil")
print("\nAn example of this table is given in example/data/dataset_metadata.csv")
print("dataset_id in the metadata table must match with the dataset id in simka distance matrices")
print("\nUse the following options to take it into consideration:")
print("\t-metadata-in: filename to a metadata table")
print("\t-metadata-variable: the name of the variable that you want to display in figures (the name of the column), for instance VARIABLE_NAME_1 in example above")
print("\nYou can now visualize simka results as heatmap, hierarchical clustering and pca in directory: " + output_dir_with_color)
print("\tcommand used:")
print("\t" + command)
