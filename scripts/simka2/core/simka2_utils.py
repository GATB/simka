
import os, time, sys, multiprocessing



class Simka2ResourceAllocator():

    MIN_MEMORY_PER_JOB = 500

    def __init__(self, isHPC, nbCores, maxMemory, maxJobs):
        self.isHPC = isHPC
        self.nbCores = nbCores
        self.maxMemory = maxMemory
        self.maxJobs = maxJobs
        #self.nbSamples = nbSamples
        #---
        self.nbSamples = None
        self.maxJobCount = None
        self.maxJobMerge = None
        self.memoryPerJob = None
        self.coresPerJob = None

        if self.nbCores == 0:
            self.nbCores = multiprocessing.cpu_count()


    def execute(self):
        if self.isHPC:
            self.execute_HPC()
        else:
            self.execute_singleNode()

    def execute_singleNode(self):

        maxjob_byCore = self.nbCores/2

        if not (self.nbSamples is None):
            maxjob_byCore = min(maxjob_byCore, self.nbSamples)

        maxjob_byCore = max(maxjob_byCore, 1)

        #maxjob_byCore = min(self.nbCores/2, self.nbSamples)

        maxjob_byMemory = self.maxMemory/Simka2ResourceAllocator.MIN_MEMORY_PER_JOB
        maxjob_byMemory = max(maxjob_byMemory, 1)

        maxJobs = min(maxjob_byCore, maxjob_byMemory)
        self.maxJobCount = maxJobs

        if self.maxJobMerge is None:
            self.maxJobMerge = self.nbCores


        self.maxJobCount = min(self.maxJobCount, self.nbCores)
        self.maxJobMerge = min(self.maxJobMerge, self.nbCores)

        self.coresPerJob = self.nbCores / self.maxJobCount
        self.coresPerJob = max(1, self.coresPerJob)

        self.memoryPerJob = self.maxMemory / self.maxJobCount
        self.memoryPerJob = max(self.memoryPerJob, Simka2ResourceAllocator.MIN_MEMORY_PER_JOB)

        self.coresPerMergeJob = self.nbCores / self.maxJobMerge
        self.coresPerMergeJob = max(1, self.coresPerMergeJob)

        #print self.coresPerJob, self.memoryPerJob, self.maxJobCount, self.maxJobMerge

    def execute_HPC(self):

        if self.maxJobs == 0:
            maxjob_byCore = self.nbCores/2
            #maxjob_byCore = min(self.nbCores/2, self.nbSamples)

            maxjob_byCore = max(maxjob_byCore, 1)

            maxjob_byMemory = self.maxMemory/Simka2ResourceAllocator.MIN_MEMORY_PER_JOB
            maxjob_byMemory = max(maxjob_byMemory, 1)

            maxJobs = min(maxjob_byCore, maxjob_byMemory)
            self.maxJobCount = maxJobs
        else:
            self.maxJobCount = self.maxJobs
            self.maxJobMerge = self.maxJobs


        self.maxJobCount = min(self.maxJobCount, self.nbCores)
        self.maxJobMerge = min(self.maxJobMerge, self.nbCores)

        self.coresPerJob = self.nbCores / self.maxJobCount
        self.coresPerJob = max(1, self.coresPerJob)

        self.memoryPerJob = self.maxMemory / self.maxJobCount
        self.memoryPerJob = max(self.memoryPerJob, Simka2ResourceAllocator.MIN_MEMORY_PER_JOB)

        self.coresPerMergeJob = self.nbCores / self.maxJobMerge
        self.coresPerMergeJob = max(1, self.coresPerMergeJob)





class ProgressBar():

    def __init__(self, max):
        self.max = max
        self.progress = 0

    def step(self, value):
        self.progress += value
        self.display()

    def display(self):
        progress_percent = float(self.progress) / float(self.max) * 100
        #---
        sys.stdout.write('\r')
        sys.stdout.write(str(round(progress_percent, 1)) + "%")

        if self.progress == self.max:
            sys.stdout.write("\n")

        sys.stdout.flush()





class JobScheduler():

    def __init__(self, maxJobs, progressTotalJobs):
        self.maxJobs = maxJobs
        self.nbJobs = 0
        self.jobQueue = []
        self.jobQueueToRemove = []

        self.progressBar = ProgressBar(progressTotalJobs)

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
                    self.progressBar.step(1)

            if isJobAvailbale:
                for checkPointFilenameToRemove in self.jobQueueToRemove:
                    self.jobQueue.remove(checkPointFilenameToRemove)

                self.jobQueueToRemove = []
                break
            else:
                time.sleep(0.5)

        #print "job finished:  ", self.nbJobs, self.maxJobs

    def join(self):
        while(self.nbJobs > 0):
            self.wait()