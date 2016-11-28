

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])



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
#os.system(command + " > /dev/null 2>&1  ")
os.system(command)
print("\nYou can now visualize simka results as heatmap, hierarchical clustering and pca in " + distance_matrices_dir)
print("\ncommand used:")
print("\t" + command)
