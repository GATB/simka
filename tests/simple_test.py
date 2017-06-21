
import sys, os, shutil, glob, gzip
os.chdir(os.path.split(os.path.realpath(__file__))[0])

suffix = " > /dev/null 2>&1"
dir = "__results__"

def clear():
	if os.path.exists("temp_output"):
		shutil.rmtree("temp_output")
	if os.path.exists("__results__"):
		shutil.rmtree("__results__")
	os.mkdir(dir)


def decompress_simka_results(dir):
	result_filenames = glob.glob(os.path.join(dir, '*.csv.gz'))
	for filename_gz in result_filenames:
		#filename_gz = result_dir + "/" + filename
		with gzip.open(filename_gz, 'rb') as f:
			outFile = open(filename_gz[:-3], "w")
			outFile.write(f.read())
			outFile.close()
			os.remove(filename_gz)

def __test_matrices(simka_vs_truth, result_dir, truth_dir):

	ok = True

	decompress_simka_results(result_dir)
	result_filenames = glob.glob(os.path.join(result_dir, '*.csv'))
	if len(result_filenames) == 0:
		print("Error: no results")
		exit(1)

	if simka_vs_truth:
		truth_filenames = glob.glob(os.path.join(truth_dir, '*.csv'))
	else: #simka vs simka
		#if result_dir+"/mat_abundance_jaccard.csv" in truth_filenames: #comparing simka results vs simka results
		#truth_filenames.remove(result_dir+"/mat_abundance_jaccard.csv") #This distance is computed from Bray Curtis distance
		decompress_simka_results(truth_dir)
		truth_filenames = glob.glob(os.path.join(truth_dir, '*.csv'))

	truth_filenames.sort()
	result_filenames.sort()

	for result_filename in result_filenames:
		distanceName = os.path.split(result_filename)[1]
		for truth_filename in truth_filenames:
			distanceName2 = os.path.split(truth_filename)[1]
			if distanceName != distanceName2: continue

			res_file = open(result_filename, "r")
			truth_file = open(truth_filename, "r")

			#print res_file, truth_file
			res_str = res_file.read()
			truth_str = truth_file.read()

			res_file.close()
			truth_file.close()

			if(res_str != truth_str):
				print("\t- TEST ERROR:    " + distanceName)
				ok = False

	return ok


def test_dists(dir):
	if(__test_matrices(True, "__results__/" + dir, "truth/" + dir)):
		print("\tOK")
	else:
		print("\tFAILED")
		sys.exit(1)


def test_parallelization():
	if(__test_matrices(False, "__results__/results_resources1", "__results__/results_resources2")):
		print("\tOK")
	else:
		print("\tFAILED")
		sys.exit(1)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------


#test k=31 t=0
clear()
print("TESTING k=31 t=0")
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_k31_t0 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 31 -abundance-min 0 -verbose 0"
print(command)
os.system(command + suffix)
test_dists("results_k31_t0")

#test k=21 t=0
clear()
print("TESTING k=21 t=0")
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_k21_t0 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -verbose 0"
print(command)
os.system(command + suffix)
test_dists("results_k21_t0")

#test k=31 t=2
clear()
print("TESTING k=31 t=2")
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_k31_t2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 31 -abundance-min 2 -verbose 0"
print(command)
os.system(command + suffix)
test_dists("results_k31_t2")

#test k=21 t=2
clear()
print("TESTING k=21 t=2")
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_k21_t2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 2 -verbose 0"
print(command)
os.system(command + suffix)
test_dists("results_k21_t2")

#test resources 1
clear()
print("TESTING parallelization")
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_resources1 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -nb-cores 20 -max-memory 4000  -verbose 0"
os.system(command + suffix)
command = "../build/bin/simka -in ../example/simka_input.txt -out ./__results__/results_resources2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -nb-cores 2 -max-memory 2000  -verbose 0"
os.system(command + suffix)
test_parallelization()

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
clear()
