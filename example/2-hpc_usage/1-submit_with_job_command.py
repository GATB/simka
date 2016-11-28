
import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])


#-----------------------------------------------------------------------------
# Set required simka parameters
#-----------------------------------------------------------------------------
simka_script_dir = "../../scripts/simka2"
input_filename = "../data/simka_input.txt"
output_dir = "simka_results"
temp_dir = "simka_temp"


#-----------------------------------------------------------------------------
# Create the command to execute
#-----------------------------------------------------------------------------
command = "python "
command += os.path.join(simka_script_dir, "simka-hpc.py")
command += " -in " + input_filename
command += " -out " + output_dir
command += " -out-tmp " + temp_dir

#-----------------------------------------------------------------------------
# Add the HPC parameters
#-----------------------------------------------------------------------------
command += " -nb-cores " + "16"
command += " -max-memory " + "8000"
command += " -max-jobs " + "10"
command += " -submit-command " + "\"qsub -pe make 8 -M 8000\""

#-----------------------------------------------------------------------------
# Run Simka
#-----------------------------------------------------------------------------
#os.system(command + " > /dev/null 2>&1  ")

#-----------------------------------------------------------------------------
# Clear temporary computation dir
#-----------------------------------------------------------------------------
if os.path.exists(temp_dir):
	shutil.rmtree(temp_dir)

print("\nThis script can only be run if your system is equiped of a job scheduling system")
print("\nHigh perfomance computing (HPC) can be achieved by using the script " + os.path.join(simka_script_dir, "simka-hpc.py"))
print("This special version of Simka provides new options specific to HPC:")
print("\t-nb-cores X: simka will use X cores PER job")
print("\t-max-memory X: simka will use a maximum of X MB of memory PER job")
print("\t-max-jobs X: simka will submit a maximum of X jobs simultaneously")
print("\t-submit-command X:")
print("\t\tthe prefix to add to commands to be able to submit them to your job scheduling system.")
print("\t\tThis is usually an executable containg the string \'sub\'")
print("\t\tAlso indicate allowed number of cores and memory per job in this command if it's possible (see example below)")

print("\n\tcommand used:")
print("\t" + command)
