
import sys, os, shutil
os.chdir(os.path.split(os.path.realpath(__file__))[0])

suffix = " > /dev/null 2>&1"
dir = "__results__"

def clear():
	if os.path.exists("temp_output"):
		shutil.rmtree("temp_output")
	if os.path.exists("__results__"):
		shutil.rmtree("__results__")
	os.mkdir(dir)


def __test_matrices(result_dir, truth_dir):

	ok = True

	result_filenames = os.listdir(result_dir)
	result_filenames.remove("mat_abundance_jaccard.csv") #This distance is computed from Bray Curtis distance
	truth_filenames = os.listdir(truth_dir)
	if "mat_abundance_jaccard.csv" in truth_filenames:
		truth_filenames.remove("mat_abundance_jaccard.csv") #This distance is computed from Bray Curtis distance

	for i in range(0, len(result_filenames)):

		res_file = open(result_dir + "/" + result_filenames[i], "r")
		truth_file = open(truth_dir + "/" + truth_filenames[i], "r")
		
		res_str = res_file.read()
		truth_str = truth_file.read()

		res_file.close()
		truth_file.close()

		if(res_str != truth_str):
			print("\t- TEST ERROR:    " + result_dir + "    " + result_filenames[i])
			ok = False

	return ok


def test_dists(dir):
	if(__test_matrices("__results__/" + dir, "truth/" + dir)):
		print("\tOK")
	else:
		print("\tFAILED")
		sys.exit(1)


def test_parallelization():
	if(__test_matrices("__results__/results_resources1", "__results__/results_resources2")):
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
