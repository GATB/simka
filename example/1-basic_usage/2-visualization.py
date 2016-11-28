

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

if not os.path.exists(distance_matrices_dir):
    print("Please run first example before: 1-simple_test.py")
    exit(1)

#-----------------------------------------------------------------------------
# Create the command to execute
#-----------------------------------------------------------------------------
command = "python "
command += visualization_dir + "/run-visualization.py"
command += " " + distance_matrices_dir

#-----------------------------------------------------------------------------
# Run visualization script
#-----------------------------------------------------------------------------
os.system(command + " > /dev/null 2>&1  ")



print("\nYou can now visualize simka results as heatmap, hierarchical clustering and pca in " + distance_matrices_dir)
print("\ncommand used:")
print("\t" + command)
