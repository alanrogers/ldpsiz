#!/usr/bin/python

###
#@file vcf2gtp.py
#@anchor vcf2gtp
#@brief Convert .vcf files into .gtp format.
#
#vcf2gtp.py, a program that converts .vcf files into .gtp format
#===============================================================
#
#    usage: vcf2gtp [options filename.vcf geneticMapFile]
#    where options may include:
#    
#      -r <x> or --recombPerBP <x>     Set recombination rate per nucleotide.
#    
#    If two filenames are provided, the first must be a vcf file and the
#    second a file of recombination rates. In this case, the -r option is
#    not allowed.
#    
#    If one filename is provided, it must be a vcf file, and the -r option
#    is mandatory.
#    
#    If no filename is provided, the program reads vcf data from standard
#    input and the -r option is mandatory.
#    
#    The program writes to standard output.
#
#If no file is provided for recombination rates, the rate must be
#specified using the `-r` or `--recombPerBP` option, and the program
#generates map distances assuming that the recombination rate between
#two sites is proportional to the number of nucleotides that separate
#them.
#
#The genetic map file must have a format like this:
#
#    position COMBINED_rate(cM/Mb) Genetic_Map(cM)
#    14431347 9.6640973708 0
#    14432618 9.7078062447 0.0122830678
#    14433624 9.7138922111 0.0220491208
#    14433659 9.7163435060 0.0223891071
#    14433758 9.7078087850 0.0233510251
#
#The first line is a list of column labels. After that, each row
#corresponds to a SNP. The first column gives the position, in
#nucleotides from the end of the chromosome. The second column gives
#the recombination rate, in cM per Mb, between that SNP and the next
#one, and the third column gives the map distance in cM. This is the
#format used by the 1000-Genomes Project
#
# @author Ryan Bohlender and Alan Rogers
# @copyright Copyright (c) 2014, Ryan Bohlender and Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".

import sys
import datetime

# Print usage message and abort
def usage(msg1):
    msg1 += "\n"
    sys.stderr.write(msg1)
    msg = \
        """
usage: vcf2gtp [options] filename.vcf [geneticMapFile]
where options may include:

  -r <x> or --recombPerBP <x>     Set recombination rate per nucleotide.

If two filenames are provided, the first must be a vcf file and the
second a file of recombination rates. In this case, the -r option is
not allowed.

If one filename is provided, it must be a vcf file, and the -r option
is mandatory.

If no filename is provided, the program reads vcf data from standard
input and the -r option is mandatory.

The program writes to standard output.
"""
    sys.stderr.write(msg)
    exit(1)

### Represents the genetic map. Map data can be provided either via a
### fixed recombination rate per Mb, or via a character, string
### representing a file name, which should be should be in the format
### provided by the 1000 Genomes Project. A Map is created like this:
###
###    mf = Map(cPerBasePair, mapFileName)
###
### Exactly one of these arguments must equal None.  If
### mapFileName!=None, the class maintains an eof variable (either
### True or False) and the following data about two adjacent map
### positions: nucpos (nucleotide position, an integer representing
### the position of the SNP on the chromosome in base-pairs), mappos
### (genetic map position in cM, a float), and cPerBasePair (local
### recombination rate per base pair).
###
### If cPerBasePair!=None, the the class retains only the value of the
### argument and an eof variable, which is always equal to True.
class Map:
    def __init__(self, cPerBasePair, mapFileName, extrapolate=False):
        self.eof = False
        self.cPerBasePair = cPerBasePair
        self.extrapolate = extrapolate
        if cPerBasePair != None and mapFileName != None:
            print "Error: Map(non-None, non-None)"
            sys.exit(1)
            return

        if cPerBasePair==None and mapFileName==None:
            print "Error: Map(None, None)"
            sys.exit(1)
            return

        if mapFileName==None:
            return

        # Set up map file
        self.mfile = open(mapFileName)

        # read header of map file and discard; set instance variables
        line = self.mfile.readline().strip().split()
        if len(line) != 3 or line[0] != "position":
            sys.stderr.write("ERR: map file has wrong format")
            sys.exit(1)
        self.nucposB=0
        self.mapposB=self.cMperMbB=0.0

        # make A and B refer to 1st two map locations
        self.next()
        self.next()

    def close(self):
        if self.cPerBasePair==None:
            self.mfile.close()
        return

    def next(self):
        """
        Move MapFile forward to next mapped SNP. If there is no more
        data, self.eof is set to True, and nothing else changes. Otherwise,
        self.eof is set to False, the values of nucposA, mapposA, and cMperMbA
        are set to the old values of nucposB, mapposB, and cMperMbB,
        and the latter three values are set from the next line of the file.
        """

        # If self.eof, or if self.cPerBasePair!=None, then this
        # function is a noop, so return.
        if self.cPerBasePair!=None or self.eof==True:
            return

        # Read new B from file. On eof, return without changing A or B.
        line = self.mfile.readline()
        if len(line) == 0: # EOF
            self.eof = True
            return
        line = line.strip().split()

        # A becomes old B; B is from file
        self.nucposA, self.mapposA, self.cMperMbA = \
            self.nucposB, self.mapposB, self.cMperMbB
        self.nucposB, self.mapposB, self.cMperMbB = \
            int(line[0]), float(line[2]), float(line[1])
        return
        
    # return map position in centimorgans
    def getMapPos(self, nucposIn):

        # If self.cPerBasePair is set, then map position in cM is a
        # linear function of nucleotide position.
        if self.cPerBasePair:
            assert mapFile==None
            return 100*nucposIn*self.cPerBasePair

        # For cases in which SNP is before 1st entry in genetic map
        if nucposIn < self.nucposA:
            if self.extrapolate:
                sys.stderr.write("extrapolating left\n")
                mBdiff = (self.nucposA - nucposIn)*1e-6
                return self.mapposA - mBdiff*self.cMperMbA
            else:
                return None

        if nucposIn == self.nucposA:
            return self.mapposA

        assert nucposIn > self.nucposA

        while self.nucposB < nucposIn and self.eof==False:
            self.next()

        if self.nucposB == nucposIn:
            return self.mapposB
        if nucposIn < self.nucposB:
            w = float(nucposIn - self.nucposA)/(self.nucposB - self.nucposA)

            # linear interpolation
            return w*self.mapposB + (1-w)*self.mapposA

        # otherwise, nucposIn is larger than last value in mapFile
        if self.extrapolate:
            sys.stderr.write("extrapolating right\n")
            assert self.nucposB < nucposIn
            mBdiff = (nucposIn - self.nucposB)*1e-6
            rtnval = self.mapposB + mBdiff*self.cMperMbB
        else:
            rtnval = None

        return rtnval

def VCFparse(line):
    """
    Utility function that takes a snp line from a vcf file
    and returns (chrom, nucpos, alleles, genotypes)
    """
    l = line.strip('\n').split('\t')
    gt = l[9:] #Skipping the mandatory fields in genotype snp data
    chrom = l[0].strip()
    nucpos = int(l[1])
    alleles = ''.join(l[3:5])
    if len(alleles) > 2:
        return (chrom,nucpos,alleles,"")  # skip loci with multiple alleles
    gtypes = []
    for ind in gt:
        ind = ind[0:3]  #the genotype and phasing info is in the first 3 positions
        assert ind[0] in "01"
        assert ind[1] in "|/"
        assert ind[2] in "01"
        if ind[1] == '|':  #Phased data. Use as is, removing pipe.
            gtypes.append(''.join(ind.split('|')))
        else:   #if unphased
            if ind[0] != ind[2]:   #if heterozygous
                gtypes.append('h')  #h replaces unphased heterozygotes
            else:
                gtypes.append(''.join(ind.split('/'))) #treat homozygotes as is
    gtypes = "".join(gtypes)
    return (chrom, nucpos, alleles, gtypes)

# BEGIN MAIN

# initialization
mapFname = None
map = None
vcfFile = None
recomb = 0.0
extrapolate = False

# Loop over command line arguments, ignoring the 0th.
# (The 0th is just the name of the program.)
i = 1
while(True):
    if i >= len(sys.argv):
        break
    elif sys.argv[i]=="-r" or sys.argv[i]=="--recombPerBP":
        i += 1
        if i >= len(sys.argv):
            usage("Missing arg to -r or --recombPerBP")
        recomb = float(sys.argv[i])
    elif sys.argv[i]=="-e" or sys.argv[i]=="--extrapolate":
        extrapolate = True
    elif sys.argv[i][0] == "-":
        usage("Unknown argument: %s" % sys.argv[i])
    else:
        # open an input file
        if vcfFile == None:
            vcfFname = sys.argv[i]
            vcfFile = open(vcfFname, "r")
        else:
            if mapFname != None:
                usage("Only 2 input files are allowed")
            mapFname = sys.argv[i]
    i += 1

# Construct Map    
if mapFname!=None:
    map = Map(None, mapFname, extrapolate)
else:
    map = Map(recomb, None, extrapolate)

if vcfFile == None:    
    vcfFile = sys.stdin
    vcfFname = "<sys.stdin>"

print "# %-19s = %s" % ("vcf2gtp.py run at", datetime.datetime.now())
print "# %-19s = %s" % ("vcf input file", vcfFname)
if mapFname!=None:
    if recomb > 0.0:
        usage("Argument -r (or --recombPerBP) can't be used if map file is provided")
    print "# %-19s = %s" % ("recombination file", mapFname)
elif  recomb > 0.0:
    print "# %-19s = %lf" % ("recombination rate", recomb)
else:
    usage("No recombination data")
    
print "# %-19s = %d" % ("ploidy", 2)

if vcfFile is sys.stdin:
    sys.stderr.write("Reading vcf data from standard input.\n")

# Read vcf header
while True:
    line = vcfFile.readline()
    if line.startswith('##ref'): #ensure reference is included
        #needed so that chr pos pairs are unique
        print "# %-19s = %s" % ("reference", line[2:].rstrip())
    elif line.startswith('#'):
        continue
    elif ',' in line:
        continue
    else:
        break

snpid = 0

# Process 1st line of vcf input. 
if line:
    chrom0, nucpos0, alleles, genotypes = VCFparse(line)
    print "# %-19s = %s" % ("chromosome", chrom0)
    print "# %-19s = %d" %("Haploid sample size",
                           len(genotypes) + genotypes.count("h"))
    print "#%9s %10s %16s %7s %s" \
        % ("snp_id", "nucpos", "mappos", "alleles", "genotypes")
    mappos = map.getMapPos(nucpos0)
    if mappos != None:
        if len(alleles) == 2:
            snpid += 1
            print "%10d %10ld %16.12lf %7s %s" % (snpid, nucpos0, mappos,
                                                  alleles, genotypes)

# Read SNP data
for line in vcfFile:
    chrom, nucpos, alleles, genotypes = VCFparse(line)
    if chrom != chrom0:
        msg = "ERR: vcf file contains multiple chromosomes: %s, %s.\n" % \
            (chrom0, chrom)
        msg += "  Detected at nucleotide position %d." % nucpos
        sys.stderr.write(msg)
        exit(1)
    if nucpos < nucpos0:
        msg = "ERR: nucleotide positions are unsorted in vcf file.\n"
        msg += "  Detected at nucleotide position %d." % nucpos
        sys.stderr.write(msg)
        exit(1)
    if nucpos == nucpos0:
        sys.stderr.write("WARNING: duplicate nucleotide positions in vcf file")
        sys.stderr.write(" Detected at nucleotide position %d"%nucpos)
    mappos = map.getMapPos(nucpos)
    if mappos == None:
        if map.eof:
            break
        else:
            continue

    nucpos0 = nucpos


    if len(alleles) == 2 and mappos >= 0.0:
        snpid += 1
        print "%10d %10ld %16.12lf %7s %s" % (snpid, nucpos, mappos,
                                              alleles, genotypes)
        
vcfFile.close()
map.close()
