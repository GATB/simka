
import os, sys
from simka2_utils import JobScheduler

#SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]

commandsFilename = sys.argv[1]
#print commandsFilename
#print checkpointFilename
#command = ""
#for i in range(2, len(sys.argv)):
#    command += sys.argv[i] + " "

commands = []
commandFile = open(commandsFilename, "r")
for command in commandFile:
    command = command.strip()
    if command == "": continue
    commands.append(command)


def jobEnd(data):
    pass

jobScheduler = JobScheduler(len(commands))


jobScheduler.start()

for command in commands:
    checkPointFilename, run_command =  command.split("|")
    #print checkPointFilename, "  -------------  ",run_command
    os.system(run_command)
    print run_command
    jobScheduler.submitJob((checkPointFilename, jobEnd, ()))

jobScheduler.join()
