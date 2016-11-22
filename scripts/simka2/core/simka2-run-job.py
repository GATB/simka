
import os, sys

checkpointFilename = sys.argv[1]
#print checkpointFilename
command = ""
for i in range(2, len(sys.argv)):
    command += sys.argv[i] + " "

command += " > /dev/null 2>&1  "
#print command

core_dump_error = 35584
ret = core_dump_error

while(ret == core_dump_error):
    ret = os.system(command)
    print ret

print ret
#Success
if ret == 0:
    f = open(checkpointFilename + "success", "w").close()
else:
    f = open(checkpointFilename + "unsuccess", "w").close()