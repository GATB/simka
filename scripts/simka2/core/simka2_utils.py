
import os, time, sys, multiprocessing

import datetime
import dateutil.relativedelta

class Simka2ResourceAllocator():

    MIN_MEMORY_PER_JOB = 500

    def __init__(self, isHPC, nbCores, maxMemory, maxJobs):
        self.isHPC = isHPC
        self.nbCores = nbCores
        self.maxMemory = maxMemory
        self.maxJobs = maxJobs
        #self.nbSamples = nbSamples
        #---
        #self.nbSamples = None
        #self.maxJobs = None
        #self.memoryPerJob = None
        #self.coresPerJob = None

        if self.nbCores == 0:
            self.nbCores = multiprocessing.cpu_count()

    def executeForCountJobs(self, nbSamplesToProcess):
        if self.isHPC:
            return self.execute_count_HPC(nbSamplesToProcess)
        else:
            return self.execute_count_singleNode(nbSamplesToProcess)

    def executeForDistanceJobs(self, nbPartitions):
        if self.isHPC:
            return self.execute_distance_HPC(nbPartitions)
        else:
            return self.execute_distance_singleNode(nbPartitions)

    def execute_count_singleNode(self, nbSamplesToProcess):

        maxjob_byCore = self.nbCores/2
        maxjob_byCore = min(maxjob_byCore, nbSamplesToProcess)
        maxjob_byCore = max(maxjob_byCore, 1)

        maxjob_byMemory = self.maxMemory/Simka2ResourceAllocator.MIN_MEMORY_PER_JOB
        maxjob_byMemory = max(maxjob_byMemory, 1)

        maxJobs = min(maxjob_byCore, maxjob_byMemory)
        maxJobs = min(maxJobs, self.nbCores)

        jobCores = self.nbCores / maxJobs
        jobCores = max(1, jobCores)

        jobMemory = self.maxMemory / maxJobs
        jobMemory = max(jobMemory, Simka2ResourceAllocator.MIN_MEMORY_PER_JOB)

        #print self.coresPerJob, self.memoryPerJob, self.maxJobCount, self.maxJobMerge
        return (maxJobs, jobCores, jobMemory)

    def execute_distance_singleNode(self, nbPartitions):

        maxJobs = 0

        if nbPartitions == -1:
            maxJobs = self.nbCores
        else:
            maxJobs = nbPartitions

        maxJobs = min(maxJobs, self.nbCores)

        jobCores = self.nbCores / maxJobs
        jobCores = max(1, jobCores)

        return (maxJobs, jobCores)

    def execute_count_HPC(self, nbSamplesToProcess):
        print("TODO")



    def execute_distance_HPC(self, nbPartitions):
        print("TODO")




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

                checkPointFilename = jobData[0]

                if os.path.exists(checkPointFilename):
                    self.jobQueueToRemove.append(jobData)
                    isJobAvailbale = True
                    self.nbJobs -= 1
                    jobData[1](jobData[2]) #Call job end method (jobData[1]) with args (jobData[2])
                    if self.progressBar is not None: self.progressBar.step(1)

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