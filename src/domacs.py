#!/usr/bin/python

###
# @page domacs.py
# @anchor domacs
# @author Alan R. Rogers
# @brief Run Gary Chen's program, "MACS"
#
# domacs.py: Use MACS to simulate a history of constant population size.
# =====================================================================
#
# Usage: domacs.py [-r or --run] output_file_basename
#
# Chen doesn't say how time is scaled. I'm guessing that units are
# 4N generations, because that is how recombination rate and mutation
# are scaled.
#
# @copyright Copyright (c) 2014, Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".


import os, sys
from ini import *

def usage():
    print "Usage: domacs.py [options] output_file_basename"
    print "  where options may include:"
    print "  -r or --run: run macs. By default, command is printed"
    print "               but not executed."
    sys.exit(1)

runprogram = False
obasename = None
for arg in sys.argv[1:]:
    if arg == "-r" or arg == "--run":
        runprogram = True
    elif arg[0] == '-':
        usage()
    else:
        if obasename == None:
            obasename = arg[:]
        else:
            usage()

if obasename==None:
    usage()

# Parameters defined here rather than in initialization file
retain = 100                # explained in MACS documentation
outputfile = obasename + ".macs"
errfile = obasename + ".macserr"
print "Writing SNPs to", outputfile
print "Writing stderr to", errfile

# name of initialization file
iniFile = "ldpsiz.ini"

# Read parameter definitions from initialization file
a, ph = readIni("ldpsiz.ini")

# Make convenience names for variables defined in initialization file
twoN0 = ph[0].twoN
if len(ph) > 1:
    twoN1 = ph[1].twoN

t0 = ph[0].t

twoNsmp = int(a["twoNsmp"])
basepairs = float(a["basepairs"])
mutation = float(a["mutation"])
recombination = float(a["recombination"])

# construct command string
cmd = "macs %s %s" % (twoNsmp, basepairs)
cmd += " -h %s" % retain
cmd += " -t %g" % (mutation * 2 * twoN0)
cmd += " -r %g" % (recombination * 2 * twoN0)

#  Add additional epochs of population history
time = 0.0
for i in range(len(ph)):
    time += ph[i].t
    if isfinite(time):
        cmd += " -eN %g %g" % (time/(2*twoN0), ph[i+1].twoN/twoN0)
    
cmd += " 2>%s 1>%s" % (errfile, outputfile)


if runprogram:
    print "Running:", cmd
    os.system(cmd)
else:
    print "Dry run:", cmd

