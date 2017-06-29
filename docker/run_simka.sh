#!/usr/bin/env bash
#
# A script to be used within a Docker container: it aims at starting a simka 
# task given some parameters.
#
# Use: ./run_simka.sh -c <command> -- <arguments>
#
#        <command>: MUST BE one of: simka, visu, test
#      <arguments>: remaining arguments passed in after <command> are passed
#                   to the appropriate simka program: 
#                      bin/simka
#                      scripts/visualization/run-visualization.py
#                   Please refer to these programs to review their expected arguments.
#
# Author: Patrick G. Durand, Inria, June 2017
# ========================================================================================
# Section: utility function declarations
# --------
# FUNCTION: display help message
function help(){
  printf "\n$0: a tool to invoke simka within a Docker container.\n\n"
  printf "usage: $0 -c <command> [arguments]\n\n"
  exit 1
}

# ========================================================================================
# Section: Main

# Prepare arguments for processing
while getopts hc: opt
do
  case "$opt" in
    c)  COMMAND="$OPTARG";;
    h)  help;;
    \?) help;;
  esac
done
shift `expr $OPTIND - 1`

# remaining arguments, if any, are supposed to be the [file ...] part of the command-line
ALL_ARGS=$@

#execute command
case "$COMMAND" in
  test)
    cd /opt/simka/example 
    ./simple_test.sh
    ;;
  simka)
    /opt/simka/bin/simka $ALL_ARGS
    ;;
  visu)
    python2.7 /opt/simka/scripts/visualization/run-visualization.py $ALL_ARGS
    ;;
esac

exit 0

