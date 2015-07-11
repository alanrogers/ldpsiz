# Where the executable files will be copied
destination := $(HOME)/bin

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions         # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include
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
lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm -lmpfr -lgmp

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

mpfr : mpfr.o
	$(CC) $(CFLAGS) -o $@ mpfr.o $(lib)
