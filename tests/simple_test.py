
import os, shutil

suffix = " &> /dev/null"
dir = "__results__"

if os.path.exists("__results__"):
	shutil.rmtree(dir)
os.mkdir(dir)

def __test_matrices(result_dir, truth_dir):
	
	ok = True
	
	result_filenames = os.listdir(result_dir)
	truth_filenames = os.listdir(truth_dir)
	
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


def test_parallelization():
	if(__test_matrices("__results__/results_resources1", "__results__/results_resources2")):
		print("\tOK")
	else:
		print("\tFAILED")


#test k=31 t=0
print("TESTING k=31 t=0")
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_k31_t0 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 31 -abundance-min 0 -verbose 0"
os.system(command + suffix)
test_dists("results_k31_t0")

#test k=21 t=0
print("TESTING k=21 t=0")
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_k21_t0 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -verbose 0"
os.system(command + suffix)
test_dists("results_k21_t0")

#test k=31 t=2
print("TESTING k=31 t=2")
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_k31_t2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 31 -abundance-min 2 -verbose 0"
os.system(command + suffix)
test_dists("results_k31_t2")

#test k=21 t=2
print("TESTING k=21 t=2")
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_k21_t2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 2 -verbose 0"
os.system(command + suffix)
test_dists("results_k21_t2")

#test resources 1
print("TESTING parallelization")
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_resources1 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -nb-cores 20 -max-memory 4000  -verbose 0"
os.system(command + suffix)
command = "../build/simka -in ../example/simka_input.txt -out ./__results__/results_resources2 -out-tmp ./temp_output -simple-dist -complex-dist -kmer-size 21 -abundance-min 0 -nb-cores 2 -max-memory 2000  -verbose 0"
os.system(command + suffix)
test_parallelization()





shutil.rmtree(dir)











