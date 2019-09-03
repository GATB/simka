#*****************************************************************************
#   SimkaMin: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
#   A tool from the GATB (Genome Assembly Tool Box)
#   Copyright (C) 2019  INRIA
#   Authors: G.Benoit, C.Lemaitre, P.Peterlongo
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************


import argparse, struct, time, datetime, sys, os, subprocess

def is_executable(bin):
    try:
        subprocess.call([bin, "-h"],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    except OSError as e:
        return(0)
    return(1)


#-------------------------------------------------------------------------------------------------------------
# ProgressBar
#-------------------------------------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------------------------------------
# ArgumentFormatterSimka
#-------------------------------------------------------------------------------------------------------------
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



#-------------------------------------------------------------------------------------------------------------
# Sketch reader
#-------------------------------------------------------------------------------------------------------------
def read_sketch_header(sketchFilename):
    f = open(sketchFilename, mode='rb')
    kmerSize = struct.unpack("B", f.read(1))[0] #B = unsigned char
    sketchSize = struct.unpack("I", f.read(4))[0] #I = unsigned int
    seed = struct.unpack("I", f.read(4))[0] #I = unsigned int
    nbDatasets = struct.unpack("I", f.read(4))[0] #I = unsigned int
    f.close()

    #u_int8_t kmerSize_;
    #file.read((char*)(&kmerSize_), sizeof(kmerSize_));
    #u_int32_t sketchSize_;
    #file.read((char*)(&sketchSize_), sizeof(sketchSize_));
    #u_int32_t seed_;
    #file.read((char*)(&seed_), sizeof(seed_));
    #u_int32_t nbDatasets_;
    #file.read((char*)(&nbDatasets_), sizeof(nbDatasets_));

    return {"kmerSize": kmerSize, "sketchSize": sketchSize, "seed": seed, "nbDatasets": nbDatasets}
