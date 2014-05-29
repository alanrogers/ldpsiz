###
# @file ini.py 
# @page ini
# @author Alan R. Rogers
# @brief Functions for objects of class Ini, which reads parameters
# from an initialization file 
#
# @copyright Copyright (c) 2014, Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".


#      Parameters         Used in
# 
#      PopHist            fitld predld
#      blocksize          eld
#      bootfilename       eld fitld
#      bootreps           eld
#      c_per_nucleotide   eld fitld
#      confidence         fitld
#      doEquilibria       predld
#      hiC               predld
#      loC               predld
#      methods            fitld predld
#      nbins              eld predld
#      nthreads           eld fitld
#      samplingInterval  eld
#      twoNsmp            fitld predld
#      u                  fitld predld
#      verbose            eld
#      windowsize_cm      eld
#
# @param[in] ifname A string, the name of the input file. For the
# format of this file, see ldpsiz.ini.
#
# @returns A dictionary containing the names and values of all
# assigned variables.
import sys

def isfinite(x):
    return (x < 1+x)

###  Represents an epoch of population history, within which the
###  population's characteristics do not change.
class Epoch:
    def __init__(self, line):
        """
        Input should be a list with two entries, each of which is a string.
        The first of these it interpreted as the length, t,
        of the Epoch in generations. The second is the number, twoN,
        of gene copies within the population during that Epoch.
        """
        if(len(line) != 2):
            print "Epoch: bad input: len(line)=%d" % len(line)
            sys.exit(1)

        ### length of epoch in generations
        self.t = float(line[0]) 

        ### haploid population size during epoch
        self.twoN = float(line[1])

    def __str__(self):
        s = "[t=%g, twoN=%g]" % (self.t, self.twoN)
        return s

def readIni(ifname):

    inPopHist = False
    ifp = open(ifname, "r")

    ph = []
    assignments = {}

    for lineno, line in enumerate(ifp):

        # strip comment
        i = line.find('#')
        if i >= 0:
            line = line[0:i]

        # strip blank lines
        line = line.strip()
        if len(line) == 0:
            continue

        # remove tabs and convert to lower case
        line = line.replace('\t',' ')
    
        if inPopHist:
            line = line.split()
            if len(line) != 2:
                print "ERR@%s:%d:PopHist lines must contain exactly two entries" % \
                    (ifname, lineno+1)
                print "line:", line
                print "len(line):", len(line)
                sys.exit(1)
            ph.append(Epoch(line))
        elif line.find('=') >= 0:  # assignment statement
            line = line.split('=')
            if len(line) != 2:
                print "Broken assignment @%s:%d" % (ifname, lineno+1)
                print "line:", line
                sys.exit(1)
            line[0] = line[0].strip()
            line[1] = line[1].strip()
            assignments[line[0]] = line[1]
        else:
            line = line.split(' ')
            line[0] = line[0].strip()
            if len(line) != 1:
                print "Broken command @%s:%d. inPopHist=%d" % (ifname, lineno+1, inPopHist)
                sys.exit(1)
            if line[0] == "PopHist":
                inPopHist = True
            else:
                assignments[line[0]] = "1"
    return assignments, ph

if __name__ == '__main__':
    a, ph = readIni("ldpsiz.ini")

    print "basepairs:", a["basepairs"]
    print "blocksize:", a["blocksize"]
    print "recombination:", a["recombination"]
    print "bootfilename:", a["bootfilename"]
    print "bootreps:", a["bootreps"]
    print "mutation:", a["mutation"]
    print "confidence:", a["confidence"]
    print "loCm:", a["loCm"]
    print "hiCm:", a["hiCm"]
    print "hiCm:", a["hiCm"]
    print "methods:", a["methods"]
    print "nbins:", a["nbins"]
    print "nthreads:", a["nthreads"]
    print "samplingInterval:", a["samplingInterval"]
    print "twoNsmp:", a["twoNsmp"]
    print "windowCm:", a["windowCm"]
    print "PopHist:"
    for i in range(len(ph)):
        s = "Epoch %2d:" % i
        print  s, ph[i]


# Local Variables:
# mode: python
# End:
