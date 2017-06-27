# *SIMKA* and *Docker*

This document explains how you can setup and use *SIMKA* within a Docker container.

## Requirements

Of course, you need to have [Docker](https://docs.docker.com/engine/installation/) installed on your system. 

We also suppose that you are familiar with [docker build](https://docs.docker.com/engine/reference/commandline/build/) and [docker run](https://docs.docker.com/engine/reference/commandline/run/) commands.

Note: this SIMKA's *Dockerfile* was made and tested using *Docker version 17* on *Mac OSX Yosemite*. However, it should work on other releases of Docker and OS (Linux, Windows or OSX). 

# How to build and run using the command-line?

## Build the container

    docker build -f Dockerfile -t simka_machine .

## Run a Simka job with sample data

    docker run --rm -i -t -v $PWD:/tmp simka_machine -c simka -- -in /opt/simka/example/simka_input.txt -out /tmp/simka_results/ -out-tmp /tmp/simka_temp_output

You should have results in ```$PWD/simka_results``` directory when Simka job is done.

## Run Simka Visualization job with provided data

    docker run --rm -i -t -v $PWD:/tmp simka_machine -c visu -- -in /tmp/simka_results/ -out /tmp/simka_results/ -pca -heatmap -tree

You should have PNG files in ```$PWD/simka_results``` directory.

## More documentation

Please refer to the documented header of the ```Dockerfile``` located in this directory.


# How to run Simka using the GoDocker platform?

## What is GoDocker?

[GoDocker](http://www.genouest.org/godocker/) is a front-end to execute Docker containers on the [Genouest](http://www.genouest.org) bioinformatics platform. An account is required to access this service.

## How to prepare a Simka job

* Log in to the GoDocker platform [here](https://godocker.genouest.org/) using your GenOuest credentials
* Click on ```Create Job``` (top-left toolbar)
* Then fill in the new job as follows:
  * Name: ```simka``` *(adapt to your needs)*
  * Description: ```simka job``` *(adapt to your needs)*
  * Tags *(leave empty)*
  * Projects *(leave value to 'default')*
  * Container image: ```pgdurand56/simka140``` *(see ```Comment 1``` below)*
  * Command: *(see ```Comment 2``` below)*
  * CPU: ```4```
  * GPU: *(leave value to '0')*
  * RAM: ```8```
  * Mount volumes: select ```home``` and/or ```omaha```  *(see ```Comment 2``` below)*
  * Advanced options: *(do not modify)*
* Click on [Submit]

### Comment 1: the Simka Docker Image pgdurand56/simka140

In this tutorial you'll use the [pgdurand56/simka140](https://hub.docker.com/r/pgdurand56/simka140/) Docker Image: this is an official Simka 1.4.0 runtime made by Genscale team member. If you want to use your own, see below.

### Comment 2: the Simka command to use

In order to use Simka Docker Image, you'll have to know that:

* GoDocker won't use the default entrypoint defined in [Simka Dockerfile](https://github.com/GATB/simka/blob/master/docker/Dockerfile). As a consequence, you do no start Simka on GoDocker as you do on the command-line.
* GoDocker enables you to access YOUR data located either in your *home directory* or in the *Omaha* storage on Genocluster machine 

#### Start a Simka data processing Job

So, here is an example of command to use while setting up a Simka job for GoDocker:

    #!/bin/bash
    /opt/simka/bin/simka -in $GODOCKER_HOME/simka/example/simka_input.txt -out $GODOCKER_HOME/simka/example/simka_results/ -out-tmp $GODOCKER_HOME/simka/example/simka_temp_output
    
In the above short script, we suppose that the data are located in the home directory of the user (denoted by variable $GODOCKER\_HOME). Simply adapt paths to your needs. If you want to use data located in Omaha, use '/omaha-beach' instead.

In this script, please DO NOT modify path: ```/opt/simka/bin/simka```. It targets the simka binary within the Simka Docker image.

#### Start a Simka visualization Job

After running a Simka data processing job, you can prepare PNG images using:

    #!/bin/bash
    python2.7 /opt/simka/scripts/visualization/run-visualization.py -in $GODOCKER_HOME/simka/example/simka_results/ -out $GODOCKER_HOME/simka/example/simka_results/ -pca -heatmap -tree
    
In this script:

* DO NOT modify path: "python2.7 /opt/simka/scripts/visualization/run-visualization.py". It targets a simka python script within the Simka Docker container.
* adapt the use of $GODOCKER\_HOME to your needs; you can also targets data located in Omaha using '/omaha-beach'


### Making your own Simka image for GoDocker

On your local computer:

    [1] cd /tmp
        git clone https://github.com/GATB/simka.git
    [2] cd simka/docker
        docker build -f Dockerfile -t simka_machine .
    [3] docker login -u <login> -p <pswd>
          (e.g. docker login -u pgdurand56 -p xxxx)
    [4] docker tag <imgID> <login>/<image_name>
          (e.g. docker tag 2520e066828a pgdurand56/simka140)
    [5] docker push <image_name>
          (e.g. docker push pgdurand56/simka140)

Steps are as follows:

    [1] get a copy of simka project
    [2] build the Simka Docker image
    [3] login to your DockerHub account
    [4] give a name to your Simka Docker Image 
        suitable for DockerHub publication
    [5] push the image to DockerHub

Now, on GoDocker use "\<login>/\<container_name>" (e.g. pgdurand56/simka140) to access your own Simka Image.

