
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
command += " -submit-command " + "\"qsub\""
command += " -submit-files " + "./job_template/job_template_sge.bash"

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

print("\nSome HPC system can only submit jobs through a job file, indicate a template of this file with the following option of simka-hpc:")
print("\t-submit-file X: simka will use this file to submit jobs")
print("\n\tIn this template submit file, you have to fill yourself the arguments that specify the allowed number of cores and memory")
print("\tSimka will then automatically add its commands to this file and submit them using the command provided in -submit-command")
print("\tSome example of submit file are provided in \"job_template\" dir")

print("\n\tcommand used:")
print("\t" + command)
