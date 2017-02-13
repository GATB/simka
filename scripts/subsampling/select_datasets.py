
"""
subsampling_stat_filename:
    output file of Simka when run with option -subsampling-setup
    path: [-out-tmp of simka]/subsampling_info.csv 

"""

import sys, os

input_filename = sys.argv[1]
subsampling_stat_filename = sys.argv[2]
output_filename = sys.argv[3]
minNbReads = int(sys.argv[4])

if subsampling_stat_filename == output_filename: exit(1)





input_file = open(input_filename, "r")
subsampling_stat_file = open(subsampling_stat_filename, "r")
subsampling_stat_file.next() #skip header
output_file = open(output_filename, "w")





datasets = {}
for line in input_file:
    line = line.strip()
    if line == "": continue

    line = line.replace(" ", "").replace("\t", "")
    id, filenames = line.split(":")

    datasets[id] = filenames




#print subsampling_stat_filename
for line in subsampling_stat_file:
    line = line.strip()
    if line == "": continue

    id, nbReads, averageReadSize = line.split(";")
    nbReads = int(nbReads)

    if nbReads >= minNbReads:
        output_file.write(id + ": " + datasets[id] + "\n")

input_file.close()
subsampling_stat_file.close()
output_file.close()
