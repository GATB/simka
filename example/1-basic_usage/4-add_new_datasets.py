
import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])


simka_script_dir = "../../scripts/simka2"
output_dir = "simka_results"
temp_dir = "simka_temp"

# Clear dirs
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)



print("\nIt is possible to update existing run of Simka with new datasets without recomputing everything again")
print("\tThis feature is useful if you know that your project will contain new datasets in the future")
print("\tThis feature can be enabled by using the following parameter of Simka:")
print("\t-keep-tmp: simka will kept temporary computation files on the disk in the (-out-tmp) dir")
print("\tBe aware that those temp files can takes a large amount of space on your disk")
print("\tEach new run of Simka must target the same -out-tmp dir")


#-----------------------------------------------------------------------------
# Run simka like in 1-simple_test.py, but keep simka database
#-----------------------------------------------------------------------------
command = "python "
command += simka_script_dir + "/simka.py"
command += " -in " + "../data/simka_input.txt"
command += " -out " + output_dir
command += " -out-tmp " + temp_dir
command += " -keep-tmp " #Keep simka database !
os.system(command + " > /dev/null 2>&1  ")

print("\nRun simka like in 1-simple_test.py, but keep simka database")
print("\tcommand used:")
print("\t" + command)

#-----------------------------------------------------------------------------
# Run simka like in 1-simple_test.py, but with another input file of datasets
#-----------------------------------------------------------------------------
command = "python "
command += simka_script_dir + "/simka.py"
command += " -in " + "../data/simka_input_2.txt"
command += " -out " + output_dir
command += " -out-tmp " + temp_dir
command += " -keep-tmp " #Keep simka database !
os.system(command + " > /dev/null 2>&1  ")

print("\nRun simka like in 1-simple_test.py, but with another input file of datasets (F and G)")
print("\tcommand used:")
print("\t" + command)


print("\nNow, exported matrices in " + output_dir + " includes F and G datasets")