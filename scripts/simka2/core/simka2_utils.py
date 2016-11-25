
import os, time, sys, multiprocessing, argparse, math

import datetime

class SimkaSettings():

	#MAX_OPEN_FILES = 1000
	MIN_FILES_TO_START_MERGE = 1000
	MAX_OPEN_FILES_PER_MERGE = 100

class SimkaCommand():

    def create_count(self):
        pass

    def create_distance(self):
        pass

    @staticmethod
    def createHPCcommand(command, isHPC, submitCommand):
        if isHPC:
            return submitCommand + " " + command
        else:
            return command

    @staticmethod
    def addHPCargs(command, args):
        if args._isHPC:
            command += " -hpc "
            command += " -max-jobs " + args._maxJobs
            command += " -submit-command " + "\"" + args.submit_command + "\""
            if args.submit_file != None:
                command += " -submit-file " + args.submit_file
        return command

class Simka2ResourceAllocator():

    MIN_MEMORY_PER_JOB = 500

    def __init__(self, isHPC, nbCores, maxMemory, maxJobs, submitCommand, submitFile):
        self.isHPC = isHPC
        self.nbCores = nbCores
        self.maxMemory = maxMemory
        self.maxJobs = maxJobs
        self.submitCommand = submitCommand
        self.submitFile = submitFile
        #print self.isHPC, self.nbCores, self.maxMemory, self.maxJobs
        #self.nbSamples = nbSamples
        #---
        #self.nbSamples = None
        #self.maxJobs = None
        #self.memoryPerJob = None
        #self.coresPerJob = None

        if self.nbCores == 0:
            self.nbCores = multiprocessing.cpu_count()

    def executeForCountJobs(self, nbSamplesToProcess, kmerSize):
        if self.isHPC:
            return self.execute_count_HPC(nbSamplesToProcess, kmerSize)
        else:
            return self.execute_count_singleNode(nbSamplesToProcess, kmerSize)

    def executeForDistanceJobs(self, nbPartitions, nbProcessedDatasets=-1, nbDatasets=-1):
        if self.isHPC:
            return self.execute_distance_HPC(nbPartitions, nbProcessedDatasets, nbDatasets)
        else:
            return self.execute_distance_singleNode(nbPartitions, nbProcessedDatasets, nbDatasets)

    def execute_count_singleNode(self, nbSamplesToProcess, kmerSize):

        maxjob_byCore = self.nbCores/2
        maxjob_byCore = max(maxjob_byCore, 1)
        #maxjob_byCore = max(maxjob_byCore, 1)

        minMemory = 0
        if kmerSize <= 15:
            minMemory = int(float(math.pow(4, kmerSize)*8)/(1<<20))
        else:
            minMemory = Simka2ResourceAllocator.MIN_MEMORY_PER_JOB
        minMemory = max(minMemory, 1)

        if self.maxMemory < minMemory:
            print("Not enough memory, you provide (" + str(self.maxMemory) + " MB), Simka need (" + str(minMemory) + "MB)")
            exit(1)

        maxjob_byMemory = self.maxMemory/minMemory
        maxjob_byMemory = max(maxjob_byMemory, 1)

        maxJobs = min(maxjob_byCore, maxjob_byMemory)
        maxJobs = min(maxJobs, self.nbCores)
        maxJobs = min(maxJobs, nbSamplesToProcess)

        jobCores = self.nbCores / maxJobs
        jobCores = max(1, jobCores)

        jobMemory = self.maxMemory / maxJobs
        jobMemory = max(jobMemory, minMemory)

        #print self.coresPerJob, self.memoryPerJob, self.maxJobCount, self.maxJobMerge
        return (maxJobs, jobCores, jobMemory)

    #nbProcessedDatasets is the number of datasets for which distance has already been computed by Simka
    #nbDatasets is the total amount of datasets in the database
    def execute_distance_singleNode(self, nbPartitions, nbProcessedDatasets=-1, nbDatasets=-1):

        maxJobs = 0

        if nbPartitions == -1:
            maxJobs = self.nbCores
        else:
            maxJobs = nbPartitions

        maxJobs = min(maxJobs, self.nbCores)

        #jobCores = self.nbCores / maxJobs
        #jobCores = max(1, jobCores)
        jobCores = 1

        #nbProcessedDatasets == -1 if this method is called by simka2-merge which not required memory
        #print "lala", nbProcessedDatasets
        if nbProcessedDatasets == -1:
            return (maxJobs, jobCores)


        nbDistances = 2 #AB-BrayCurtis and PA-Jaccard

        N = nbDatasets
        Nold = nbProcessedDatasets
        Nnew = nbDatasets-nbProcessedDatasets
        minMemoryPerCore = 8*N*Nnew * nbDistances

        maxCore_byMemory = int((self.maxMemory*1<<20) / minMemoryPerCore)
        print "min memory per core: ", minMemoryPerCore, N, nbProcessedDatasets
        print "max cores by memory: ", maxCore_byMemory

        N = Nnew

        #Memory required by Simka2-distance is (No = nbProcessedDatasets), 8 Bytes for storing any counts
        #Mem = 8NN + 8NNo
        #So given available memory, we can resolve this quadratic equation to know how much datasets we can process N
        #8NN + 8NNo - Mem = 0
        #a=8 b=8No c=-Mem
        if maxCore_byMemory <= 0:
            availableMemory = self.maxMemory * (1<<20) / nbDistances
            a = 8
            b = 8*nbProcessedDatasets
            c = -availableMemory

            d = b**2-4*a*c

            if d < 0:
                print "Not enougth memory"
                exit(1)
            elif d == 0:
                x1 = -b / (2*a)
            else: # if d > 0
                x1 = (-b + math.sqrt(d)) / (2*a)
                x2 = (-b - math.sqrt(d)) / (2*a)

            #print x1, x2
            N = int(x1)
            #print N

            if N <= 0:
                print("Not enough memory")
                exit(1)

            maxCore_byMemory = 1

        maxJobs = min(maxJobs, maxCore_byMemory)
        maxJobs = max(maxJobs, 1)

        return (maxJobs, jobCores, N)

    def execute_count_HPC(self, nbSamplesToProcess, kmerSize):

        #maxCountJob = self.maxJobs
        minMemory = 0
        if kmerSize <= 15:
            minMemory = int(float(math.pow(4, kmerSize)*8)/(1<<20))
        else:
            minMemory = Simka2ResourceAllocator.MIN_MEMORY_PER_JOB

        jobCores = self.nbCores
        jobMemory = self.maxMemory

        if jobMemory < minMemory:
            raise Exception("Not enough memory per jobs, you provide (" + str(jobMemory) + " MB), Simka need (" + str(minMemory) + ")")

        #maxjob_byCore = self.maxJobs/2
        #maxjob_byCore = min(maxjob_byCore, nbSamplesToProcess)
        #maxjob_byCore = max(maxjob_byCore, 1)

        #maxjob_byMemory = self.maxMemory/Simka2ResourceAllocator.MIN_MEMORY_PER_JOB
        #maxjob_byMemory = max(maxjob_byMemory, 1)

        #maxJobs = min(maxjob_byCore, maxjob_byMemory)
        #maxJobs = min(maxJobs, self.nbCores)
        maxJobs = min(self.maxJobs, nbSamplesToProcess)

        #jobCores = self.nbCores / maxJobs
        #jobCores = max(1, jobCores)

        #jobMemory = self.maxMemory / maxJobs
        #jobMemory = max(jobMemory, Simka2ResourceAllocator.MIN_MEMORY_PER_JOB)

        #print self.coresPerJob, self.memoryPerJob, self.maxJobCount, self.maxJobMerge
        return (maxJobs, jobCores, jobMemory)



    def execute_distance_HPC(self, nbPartitions, nbProcessedDatasets=-1, nbDatasets=-1):
        """
        maxJobs = 0

        if nbPartitions == -1:
            maxJobs = self.maxJobs
        else:
            maxJobs = nbPartitions

        maxJobs = min(maxJobs, self.maxJobs)

        jobCores = self.nbCores

        return (maxJobs, jobCores)
        """

        maxJobs = 0

        if nbPartitions == -1:
            maxJobs = self.maxJobs
        else:
            maxJobs = nbPartitions

        maxJobs = min(maxJobs, self.maxJobs)

        #jobCores = self.nbCores / maxJobs
        #jobCores = max(1, jobCores)
        jobCores = self.nbCores

        #nbProcessedDatasets == -1 if this method is called by simka2-merge which not required memory
        #print "lala", nbProcessedDatasets
        if nbProcessedDatasets == -1:
            return (maxJobs, jobCores)


        nbDistances = 2 #AB-BrayCurtis and PA-Jaccard

        N = nbDatasets
        Nold = nbProcessedDatasets
        Nnew = nbDatasets-nbProcessedDatasets
        minMemoryPerCore = 8*N*Nnew * nbDistances

        maxCore_byMemory = int((self.maxMemory*1<<20) / minMemoryPerCore)
        print "min memory per core: ", minMemoryPerCore, N, nbProcessedDatasets
        print "max cores by memory: ", maxCore_byMemory

        N = Nnew

        #Memory required by Simka2-distance is (No = nbProcessedDatasets), 8 Bytes for storing any counts
        #Mem = 8NN + 8NNo
        #So given available memory, we can resolve this quadratic equation to know how much datasets we can process N
        #8NN + 8NNo - Mem = 0
        #a=8 b=8No c=-Mem
        if maxCore_byMemory <= 0:
            availableMemory = self.maxMemory * (1<<20) / nbDistances
            a = 8
            b = 8*nbProcessedDatasets
            c = -availableMemory

            d = b**2-4*a*c

            if d < 0:
                print "Not enougth memory"
                exit(1)
            elif d == 0:
                x1 = -b / (2*a)
            else: # if d > 0
                x1 = (-b + math.sqrt(d)) / (2*a)
                x2 = (-b - math.sqrt(d)) / (2*a)

            #print x1, x2
            N = int(x1)
            #print N

            if N <= 0:
                print("Not enough memory")
                exit(1)

            maxCore_byMemory = 1

        jobCores = min(jobCores, maxCore_byMemory)
        jobCores = max(jobCores, 1)

        return (maxJobs, jobCores, N)




class ProgressBar():

    def __init__(self, text, max):
        self.text = text
        self.max = max
        self.progress = 0
        self.start_time = 0

    def start(self):
        self.progress = 0
        self.start_time = time.time()
        self.display()

    def step(self, value):
        self.progress += value
        self.display()

    def display(self):
        progress_percent = float(self.progress) / float(self.max) * 100

        duration = int(time.time() - self.start_time)
        duration_str = str(datetime.timedelta(seconds=duration))
        #---
        sys.stdout.write('\r')
        sys.stdout.write("[" + str(round(progress_percent, 1)) + "%] " +  self.text + "    [Time: " + duration_str + "]")

        if self.progress == self.max:
            sys.stdout.write("\n")

        sys.stdout.flush()





class JobScheduler():

    def __init__(self, maxJobs, progressBar=None):
        self.maxJobs = maxJobs
        self.nbJobs = 0
        self.jobQueue = []
        self.jobQueueToRemove = []
        #---
        self.progressBar = progressBar

    def start(self):
        self.nbJobs = 0
        self.jobQueue = []
        self.jobQueueToRemove = []
        if self.progressBar is not None: self.progressBar.start()

    #jobData specif: (checkPointFilemane, endJobMethod, (endJobMethodArgs, ...))
    def submitJob(self, jobData):
        self.nbJobs += 1
        #print self.nbJobs, self.maxJobs
        self.jobQueue.append(jobData)

        #print self.nbJobs, self.maxJobs
        if self.nbJobs >= self.maxJobs:
            self.wait()



    def wait(self):

        #print("waiting")
        while True:

            #print self.jobQueue
            isJobAvailbale = False

            for jobData in self.jobQueue:

                successCheckPointFilename = jobData[0] + "success"
                unsuccessCheckPointFilename = jobData[0] + "unsuccess"

                if os.path.exists(successCheckPointFilename):
                    self.jobQueueToRemove.append(jobData)
                    isJobAvailbale = True
                    self.nbJobs -= 1
                    jobData[1](jobData[2]) #Call job end method (jobData[1]) with args (jobData[2])
                    if self.progressBar is not None: self.progressBar.step(1)
                elif os.path.exists(unsuccessCheckPointFilename):
                    print("A job failed, exiting simka")
                    exit(1)

            if isJobAvailbale:
                for checkPointFilenameToRemove in self.jobQueueToRemove:
                    self.jobQueue.remove(checkPointFilenameToRemove)

                self.jobQueueToRemove = []
                break
            else:
                time.sleep(0.1)

        #print "job finished:  ", self.nbJobs, self.maxJobs

    def join(self):
        while(self.nbJobs > 0):
            self.wait()




















class SimkaParser(argparse.ArgumentParser):

    def error(self, message):
        print("")
        sys.stderr.write('error: %s\n' % message)
        print("")
        self.print_help()
        sys.exit(2)


class ArgumentFormatterSimka(argparse.HelpFormatter):


    #def _fill_text(self, text, width, indent):
    #    return ''.join([indent + line for line in text.splitlines(True)])
    def _split_lines(self, text, width):
        return text.splitlines()

    #remove default args layout
    def _format_args(self, action, default_metavar):
        result = ""
        return result

    #Remove "usage: ..." header
    def _format_usage(self, usage, actions, groups, prefix):
        return ""


    #Changed layout of each item
    def _get_help_string(self, action):

        text = ""

        if type(action) == argparse._StoreAction:
            text =  "(1 arg) :    " + action.help
        elif type(action) == argparse._StoreTrueAction:
            text =  "(0 arg) :    " + action.help

        if type(action) == argparse._StoreAction and action.default != None:
            text += " [Default: " + str(action.default) + "]"
        #print type(action), action
        #print action
        #return "-5-"
        #return action.help
        if text != "":
            return text

        return "__none__"

    #Hack for removing useless "optional arguments:" section
    def _join_parts(self, part_strings):
        #print part_strings
        return ''.join([part
                        for part in part_strings
                        if part and part is not argparse.SUPPRESS and not "optional arguments:" in part and not "__none__" in part and not "--help" in part])


