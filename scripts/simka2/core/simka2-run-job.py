
import os, sys

checkpointFilename = sys.argv[1]
#print checkpointFilename
command = ""
for i in range(2, len(sys.argv)):
    command += sys.argv[i] + " "

#command += " > /dev/null 2>&1  "
#print command

#seg_fault_error = 11
#core_dump_error = 35584
#core_dump_error2 = 32256
#ret = core_dump_error

#while(ret == core_dump_error or ret == core_dump_error2 or ret == seg_fault_error):
print("loooool")
ret = os.system(command)
#print ret

#print ret
#Success
if ret == 0:
    f = open(checkpointFilename + "success", "w").close()
else:
    f = open(checkpointFilename + "unsuccess", "w").close()