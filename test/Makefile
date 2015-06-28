vpath %.c ../src

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src

tests := xassign xboot xchain xem xfileindex xhayes xhill \
   xini xjobqueue xmisc xmodel xpophist xreadgtp xsnp \
   xspectab xstrobeck xsums xtabulation xtokenizer xwindow \
   xmatcoal xpolya xespectrum xtfespectrum xarray

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

test : $(tests)
	-./xarray
	-./xassign
	-./xboot
	-./xchain
	-./xem
	-./xespectrum
	-./xfileindex
	-./xhayes
	-./xhill
	-./xini
	-./xjobqueue
	-./xmatcoal
	-./xmisc
	-./xmodel
	-./xpolya
	-./xpophist
	-./xreadgtp
	-./xsnp
	-./xspectab
	-./xstrobeck
	-./xsums
	-./xtabulation
	-./xtfespectrum
	-./xtokenizer 
	-./xwindow
	@echo "FINISHED UNIT TESTS."

XARRAY := xarray.o array.o misc.o sums.o
xarray : $(XARRAY)
	$(CC) $(CFLAGS) -o $@ $(XARRAY) $(lib)

XASSIGN := xassign.o assign.o misc.o tokenizer.o sums.o
xassign : $(XASSIGN)
	$(CC) $(CFLAGS) -o $@ $(XASSIGN) $(lib)

# test boot.c
XBOOT := xboot.o boot.o misc.o tabulation.o window.o snp.o readgtp.o \
	fileindex.o sums.o tokenizer.o assign.o em.o spectab.o array.o
xboot : $(XBOOT)
	$(CC) $(CFLAGS) -o $@ $(XBOOT) $(lib)

XCHAIN := xchain.o chain.o misc.o pophist.o model.o strobeck.o hayes.o \
   hill.o tokenizer.o sums.o
xchain : $(XCHAIN)
	$(CC) $(CFLAGS) -o $@ $(XCHAIN) $(lib)

XEM := xem.o em.o misc.o sums.o
xem : $(XEM)
	$(CC) $(CFLAGS) -o $@ $(XEM) $(lib)

# test fileindex.c
XFILEINDEX := xfileindex.o fileindex.o misc.o readgtp.o tokenizer.o assign.o \
    sums.o
xfileindex : $(XFILEINDEX)
	$(CC) $(CFLAGS) -o $@ $(XFILEINDEX) $(lib)

XHAYES := xhayes.o hayes.o misc.o pophist.o model.o strobeck.o tokenizer.o \
   hill.o sums.o
xhayes : $(XHAYES)
	$(CC) $(CFLAGS) -o $@ $(XHAYES) $(lib)

XHILL := xhill.o hill.o misc.o sums.o pophist.o model.o strobeck.o \
   hayes.o tokenizer.o
xhill : $(XHILL)
	$(CC) $(CFLAGS) -o $@ $(XHILL) $(lib)

XINI := xini.o ini.o misc.o sums.o tokenizer.o pophist.o assign.o
xini : $(XINI)
	$(CC) $(CFLAGS) -o $@ $(XINI) $(lib)

XJOBQUEUE := xjobqueue.o jobqueue.o misc.o sums.o
xjobqueue : $(XJOBQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XJOBQUEUE) $(lib)

XMISC := xmisc.o misc.o sums.o
xmisc : $(XMISC)
	$(CC) $(CFLAGS) -o $@ $(XMISC) $(lib)

XMODEL := xmodel.o model.o misc.o sums.o strobeck.o hill.o hayes.o \
  tokenizer.o pophist.o
xmodel : $(XMODEL)
	$(CC) $(CFLAGS) -o $@ $(XMODEL) $(lib)

XPOPHIST := xpophist.o pophist.o misc.o sums.o
xpophist : $(XPOPHIST)
	$(CC) $(CFLAGS) -o $@ $(XPOPHIST) $(lib)

XSPECTAB := xspectab.o spectab.o misc.o sums.o array.o
xspectab : $(XSPECTAB)
	$(CC) $(CFLAGS) -o $@ $(XSPECTAB) $(lib)

XREADGTP := xreadgtp.o readgtp.o misc.o sums.o tokenizer.o assign.o
xreadgtp : $(XREADGTP)
	$(CC) $(CFLAGS) -o $@ $(XREADGTP) $(lib)

XSTROBECK := xstrobeck.o strobeck.o misc.o sums.o pophist.o model.o \
  hill.o hayes.o tokenizer.o
xstrobeck : $(XSTROBECK)
	$(CC) $(CFLAGS) -o $@ $(XSTROBECK) $(lib)

XSUMS := xsums.o sums.o misc.o
xsums : $(XSUMS)
	$(CC) $(CFLAGS) -o $@ $(XSUMS) $(lib)

XTABULATION := xtabulation.o tabulation.o misc.o window.o snp.o sums.o \
    boot.o readgtp.o fileindex.o tokenizer.o assign.o em.o spectab.o \
    array.o
xtabulation : $(XTABULATION)
	$(CC) $(CFLAGS) -o $@ $(XTABULATION) $(lib)

XTOKENIZER := xtokenizer.o tokenizer.o misc.o sums.o
xtokenizer : $(XTOKENIZER)
	$(CC) $(CFLAGS) -o $@ $(XTOKENIZER) $(lib)

XWINDOW := xwindow.o window.o snp.o boot.o sums.o misc.o tabulation.o \
   readgtp.o tokenizer.o assign.o em.o spectab.o array.o 
xwindow : $(XWINDOW)
	$(CC) $(CFLAGS) -o $@ $(XWINDOW) $(lib)

XSNP := xsnp.o snp.o boot.o sums.o misc.o em.o tabulation.o \
  spectab.o array.o
xsnp : $(XSNP)
	$(CC) $(CFLAGS) -o $@ $(XSNP) $(lib)

TRY := try.o
try : $(TRY)
	$(CC) $(CFLAGS) -o $@ $(TRY) $(lib)

# test parallel reads by multiple threads
TSTPREAD := tstpread.o jobqueue.o misc.o
tstpread: $(TSTPREAD)
	$(CC) $(CFLAGS) -o $@ $(TSTPREAD) $(lib)

XMATCOAL := xmatcoal.o pophist.o matcoal.o misc.o sums.o
xmatcoal : $(XMATCOAL)
	$(CC) $(CFLAGS) -o $@ $(XMATCOAL) $(lib)

XPOLYA := xpolya.o polya.o misc.o sums.o
xpolya : $(XPOLYA)
	$(CC) $(CFLAGS) -o $@ $(XPOLYA) $(lib)

XESPECTRUM := xespectrum.o espectrum.o polya.o matcoal.o misc.o sums.o pophist.o
xespectrum : $(XESPECTRUM)
	$(CC) $(CFLAGS) -o $@ $(XESPECTRUM) $(lib)

XTFESPECTRUM := xtfespectrum.o tfespectrum.o espectrum.o polya.o matcoal.o \
  misc.o pophist.o spectab.o array.o sums.o
xtfespectrum : $(XTFESPECTRUM)
	$(CC) $(CFLAGS) -o $@ $(XTFESPECTRUM) $(lib)

install: $(alltargets)
	cp -p $(alltargets) $(destination)

# Make dependencies file
depend : *.c
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(ctargets) $(ctests) core.* vgcore.*

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

