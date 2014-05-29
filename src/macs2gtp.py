#!/usr/bin/python
###
#@page macs2gtp.py
#@anchor macs2gtp
#@file macs2gtp.py
#@brief Parse output of macs, produce gtp format, the input for eld.
#
#macs2gtp.py, a program that converts `macs` output into .gtp format 
#===================================================================
#
#
#    Usage: macs2gtp.py [options] \<input file\> 
#    where options may include:
#      -R \<x\> : set recombination rate for adjacent sites
#
# For example,
#
#     macs2gtp.py -R 1e-8 macs.out > macs.gtp
#
# Would read file `macs.out`, generate map distances assuming a
# recombination rate of 1e-8 per nucleotide, and write the output
# into macs.gtp.
#
# If the recombination rate is not specified, `macs2gtp.py` tries to get
# it from the initialization file `ldpsiz.ini`. 
#
# @copyright Copyright (c) 2014, Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".

import sys, datetime
from ini import *


### Print usage message and exit
def usage():
   print >> sys.stderr, "Usage: macs2gtp.py [options] <input file>"
   print >> sys.stderr, "where options may include:"
   print >> sys.stderr, " -R <x> : set recombination rate for adjacent sites"
   sys.exit(1)

recombination = None    

# Read parameter definitions from initialization file
a, ph = readIni("ldpsiz.ini")
twoN0 = ph[0].twoN

if "recombination" in a:    
    recombination = float(a["recombination"])

# Command line arguments
i = 1
while(True):
    if i >= len(sys.argv):
        break
    if sys.argv[i] == "-R":
        i += 1
        try:
            recombination = float(sys.argv[i])
        except:
            usage()
    elif sys.argv[i][0] == '-':
       usage()
    else:
        datafile = sys.argv[i]
    i += 1
    

if recombination == None:
    print "Error: no recombination rate"
    exit(1)

print "# %-35s =" % "data source",
for word in sys.argv:
    print "%s" % word.strip(),
print

infile = open(datafile)

# Skip comments at top of file
while(True):
    line = infile.readline()
    line = line.strip()
    if len(line)>0 and line[0] != '#':
        break

# First line of input has MACS command line
line = line.split()

if(len(line) < 3 or line[0] != "COMMAND:" or line[1] != "macs"):
    print >> sys.stderr, "1st line looks wrong:"
    for word in line:
        print >> sys.stderr, "%s" % word,
    print
    sys.exit(1)


# Find recombination rate per 4N nucleotides
try:
    i = line.index("-r")
    i += 1
except:
    i = len(line)

if i < len(line):
    fourNc = float(line[i])

    # Make sure command line agrees with ldpsiz.ini
    err = abs(fourNc - recombination*2*twoN0)/fourNc
    if err > 1e-10:
       print "MACS command not consistent with ldpsiz.ini in current directory"
       print "  According to MACS command, 4Nc=%g" % fourNc
       print "  In ldpsiz.ini, c=%g and 2N0=%g, so 4Nc=%g" \
           % (recombination, twoN0, recombination*2*twoN0)
       print "  relerr=%e" % err
       exit(1)

print "# %-35s = %s" % ("time", datetime.datetime.now())
print "# %-35s =" % "macs cmd",
for i in range(1,len(line)):
    print "%s" % line[i],
print

nnucleotides = int(round(float(line[3])))

# Second line has seed, which we will ignore
line = infile.readline().split()
if(len(line) != 2 or line[0] != "SEED:"):
    print "2nd line looks wrong:"
    for word in line:
        print "%s" % word,
    print
    sys.exit(1)

print "# %-35s = %lg" % ("recombination rate per nucleotide", recombination)
print "# %-35s = %d" % ("nnucleotides", nnucleotides)
print "# %-35s = %d" % ("ploidy", 1)

lineno = 0
for line in infile:
    line = line.strip().split()
    if(line=="" or line[0] != "SITE:"):
        continue
    if lineno==0:
        print "# %-35s = %d" % ("Haploid sample size", len(line[3]))
        print "#%9s %10s %14s %7s %s" \
            % ("snp_id", "nucpos", "mappos", "alleles", "genotypes")
    nfields = len(line)
    nucpos = int(round(nnucleotides*float(line[2])))
    centimorgan = nucpos*recombination*100.0
    print "%10ld %10ld %14.12lf %7s %s" % (lineno,
                                        nucpos,
                                        centimorgan,
                                        "01", line[nfields-1])
    lineno += 1
