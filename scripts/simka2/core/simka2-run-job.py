
import os, sys

checkpointFilename = sys.argv[1]
#print checkpointFilename
command = ""
for i in range(2, len(sys.argv)):
    command += sys.argv[i] + " "

command += " > /dev/null 2>&1  "
#print command
ret = os.system(command)
if ret != 0:
    print ret
else:
    f = open(checkpointFilename, "w").close()
