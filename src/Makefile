# Where the executable files will be copied
destination := $(HOME)/bin

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions         # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include

pytargets := macs2gtp.py vcf2gtp.py ini.py domacs.py 
ctargets := predld obsld fithayes ms2gtp boot2tbl fwdsim mcmcld sald \
   tabfboot
alltargets := $(ctargets) $(pytargets)

CC := gcc

# Removed because they generate excessive warnings
# -Wconversion \
# -Wfloat-equal \

# Flags to determine the warning messages issued by the compiler
warn := \
 -Wall \
 -Wcast-align \
 -Wcast-qual \
 -Wmissing-declarations \
 -Wmissing-prototypes \
 -Wnested-externs \
 -Wpointer-arith \
 -Wstrict-prototypes \
 -Wno-unused-parameter \
 -Wno-unused-function \
 -Wshadow \
 -Wundef \
 -Wwrite-strings

CFLAGS := -g -std=gnu99 $(warn) $(incl) $(opt) $(prof) $(osargs)
lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

all : $(ctargets)

PREDLD := predld.o pophist.o hayes.o hill.o strobeck.o misc.o \
          tokenizer.o model.o ini.o assign.o
predld : $(PREDLD)
	$(CC) $(CFLAGS) -o $@ $(PREDLD) $(lib)

OBSLD := obsld.o readgtp.o misc.o tabulation.o window.o assign.o \
         fileindex.o boot.o sums.o tokenizer.o ini.o pophist.o \
         em.o
obsld : $(OBSLD)
	$(CC) $(CFLAGS) -o $@ $(OBSLD) $(lib)

# use MCMCMC to fit pophist parameters to data
MCMCLD := mcmcld.o pophist.o hill.o strobeck.o misc.o \
    tokenizer.o tabulation.o model.o ini.o hayes.o assign.o chain.o 
mcmcld : $(MCMCLD)
	$(CC) $(CFLAGS) -o $@ $(MCMCLD) $(lib)

# fit pophist parameters to data using method of Hayes et al.
FITHAYES := fithayes.o misc.o pophist.o tokenizer.o ini.o hayes.o \
    assign.o
fithayes : $(FITHAYES)
	$(CC) $(CFLAGS) -o $@ $(FITHAYES) $(lib)

# convert output of ms into gtp format
MS2GTP := ms2gtp.o misc.o tokenizer.o ini.o assign.o pophist.o
ms2gtp : $(MS2GTP)
	$(CC) $(CFLAGS) -o $@ $(MS2GTP) $(lib)

# convert boot file into a table
BOOT2TBL := boot2tbl.o misc.o boot.o tabulation.o
boot2tbl : $(BOOT2TBL)
	$(CC) $(CFLAGS) -o $@ $(BOOT2TBL) $(lib)

# forward two-locus simulation
FWDSIM := fwdsim.o
fwdsim : $(FWDSIM)
	$(CC) $(CFLAGS) -o $@ $(FWDSIM) $(lib)

TRY := try.o
try : $(TRY)
	$(CC) $(CFLAGS) -o $@ $(TRY) $(lib)

SALD := sald.o sasimplex.o model.o misc.o hayes.o hill.o strobeck.o \
pophist.o tokenizer.o ini.o assign.o jobqueue.o boot.o tabulation.o \
rusage.o locref.o annealsched.o 
sald: $(SALD)
	$(CC) $(CFLAGS) -o $@ $(SALD) $(lib)

TABFBOOT := tabfboot.o misc.o pophist.o boot.o tokenizer.o \
  tabulation.o
tabfboot: $(TABFBOOT)
	$(CC) $(CFLAGS) -o $@ $(TABFBOOT) $(lib)

install: $(alltargets)
	cp -p $(alltargets) $(destination)

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(ctargets) core.* vgcore.*

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean
