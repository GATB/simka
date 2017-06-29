#########################################################################################
#
#  Docker file for Simka project.
#
#  It prepares a Docker container to run Simka jobs: 
#
#    - bin/simka: computing simka results from sequencing data
#    - scripts/visualization/run-visualization.py: making images from results of
#      bin/simka.
#
#########################################################################################
#
# == Docker build command:
#
#    docker build -f Dockerfile -t simka_machine .
#
# == Docker test command:
#
#    docker run --rm -i -t simka_machine -c test
#
#    -> you should see a simka test with some provided data.
#
# == Running a Simka job:
#
#    docker run --rm -i -t simka_machine -c <command> -- <args>
#
#    where:
#        <command>: MUST BE one of: simka, visu, test
#      <arguments>: remaining arguments passed in after <command> are passed
#                   to the appropriate simka program: 
#                      - simka: will run 'bin/simka' within the container
#                      - visu: will run 'scripts/visualization/run-visualization.py'
#                              within the container
#                   Please refer to these programs to review their expected arguments.
#                   See https://github.com/GATB/simka
#
# == Sample Simka job with provided data:
#    
#    docker run --rm -i -t -v $PWD:/tmp simka_machine -c simka -- -in /opt/simka/example/simka_input.txt -out /tmp/simka_results/ -out-tmp /tmp/simka_temp_output
#
#    -> you should have results in $PWD/simka_results directory when Simka job is done.
#
#    This command-line line explained:
#
#    docker run                                 [1]
#       --rm                                    [2]
#       -i -t                                   [3]
#       -v $PWD:/tmp                            [4]
#       simka_machine                           [5] 
#       -c simka                                [6]
#       --                                      [7]
#       -in /opt/simka/example/simka_input.txt  [8]
#       -out /tmp/simka_results/                [9]
#       -out-tmp /tmp/simka_temp_output         [10]
#
#       [1]-[5]: Docker arguments
#       [6]-[7]: simka container's invoker program
#       [8]-[10]: 'bin/simka' arguments
#
#       [1]: start Docker container
#       [2]: destroy container when Docker finishes
#            (it does NOT delete the 'simka_machine' image)
#       [3]: start an interactive job 
#            (for instance, you'll see messages on stdout, if any)
#       [4]: mount a volume. This is required to get the results from Simka.
#            Here, we say that current local directory will be viewed as '/tmp'
#            from the inside of the container. 
#       [5]: tell Docker which image to start: the 'simka_machine' of course.
#       [6]: ask to start the simka program. Other option is to start the 
#            'visu' task (see below). See companion file 'run_simka.sh' for
#            more information.
#       [7]: '--' is required to separate arguments [6] from the rest of the
#            command line
#       [8]: the data file to process with simka. Here we use a data file
#            provided with the simka software to test it.
#       [9]: tells simka where to put results. Of course, simka will write 
#            within /tmp directory inside the container. However, since we
#            have directive [4], data writing is actually done in $PWD, i.e.
#            a local directory.
#       [10]: tells simka where to put temporary files. 
#
# == Sample Simka Visualization job with provided data
#
#    After running the previous command, you can do this:
#
#    docker run --rm -i -t -v $PWD:/tmp simka_machine -c visu -- -in /tmp/simka_results/ -out /tmp/simka_results/ -pca -heatmap -tree
#
#    -> you should have PNG files in $PWD/simka_results directory.
#    
# == Additional notes
# 
#   Root access inside the container:
#
#     - if running: docker exec -it simka_machine bash
#
#     - if not yet running: docker run --rm -i -t simka_machine bash
#
#########################################################################################

# Simka binary available on Github (see below) is built using a 
# Debian 8 (jessie) based system on Inria Jenkins CI platform
FROM debian:jessie

# who to blame?
MAINTAINER Patrick Durand patrick.durand@inria.fr

# ###
#    We always use the latest official SIMKA release.
#
ENV SIMKA_VERSION=1.4.0
    
# ###
#     Package installation and configuration
#
RUN apt-get update && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends curl python2.7 r-base \
    && apt-get clean

# ###
#     SIMKA installation: get the binary release from Github mirror.
#
RUN cd /opt \
    && export SIMKA_TGZ=simka-v${SIMKA_VERSION}-bin-Linux.tar.gz \
    && export GIT_URL=https://github.com/GATB/simka/releases/download \
    && export SIMKA_URL=${GIT_URL}/v${SIMKA_VERSION}/${SIMKA_TGZ} \
    && curl -ksL ${SIMKA_URL} | tar xz \
    && rm -f ${SIMKA_TGZ} \
    && mv simka-v${SIMKA_VERSION}-bin-Linux simka \
    && cd simka/bin \
    && chmod +x simka* \
    && cd ../example \
    && chmod +x *.sh \
    && ./simple_test.sh

COPY run_simka.sh /opt/simka

# Fix: ensure script has exec permission
RUN chmod +x /opt/simka/run_simka.sh

# ###
#     Start simka. 
#
ENTRYPOINT ["/opt/simka/run_simka.sh"]


