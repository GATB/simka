
import os, sys

print("\n\n------\n" + str(sys.argv))

checkpointFilename = sys.argv[1]
jobType = sys.argv[2]
#print checkpointFilename
command = ""
for i in range(3, len(sys.argv)):
    command += sys.argv[i] + " "

print("\n\n------\n" + command)
#command += " > /dev/null 2>&1  "
#print command

#seg_fault_error = 11
#core_dump_error = 35584
#core_dump_error2 = 32256
#ret = core_dump_error

nbTry = 0
while(True):
    print("\n\nExec\n" + command)
    ret = os.system(command)
    print("Ret: " + str(ret))

    if ret == 0: break

    nbTry += 1
    if nbTry >= 3:
        break

    print("Nb try: " + str(nbTry))

print("End: " + str(ret))

#Success
if ret == 0 or jobType != "count":
    f = open(checkpointFilename + "success", "w").close()
else:
    f = open(checkpointFilename + "unsuccess", "w").close()