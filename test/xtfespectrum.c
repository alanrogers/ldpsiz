/**
   @file xtfespectrum.c
   @brief Test class TFESpectrum
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2015, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "tfespectrum.h"
#include "polya.h"
#include "pophist.h"
#include "misc.h"
#include "spectab.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>

int main(int argc, char **argv) {
    unsigned nSamples=10u;
    int verbose = 0, ok=1;
    unsigned truncSFS = 2;
	unsigned i, j;
	double errTol = 0.00001;
    for(i=1; i<argc; ++i) {
		if(strcmp(argv[i], "-v") == 0)
			verbose = 1;
		else if(strcmp(argv[i], "-n") == 0) {
			if(++i >= argc) {
				eprintf("%s:%d: missing arg for -n\n",
						__FILE__, __LINE__);
			}
			nSamples = strtol(argv[i], 0, 10);
		}
    }

    // Polya is expensive, so allocate once in main.
	Polya *polya = Polya_new(nSamples);
	checkmem(polya, __FILE__, __LINE__);

	// pophist has 1 epoch with pop size twoN=1.
    EpochLink  *linkedList = EpochLink_new(NULL, HUGE_VAL, 1.0);
    PopHist *ph = PopHist_fromEpochLink(linkedList);
    TFESpectrum *tfes = TFESpectrum_new(nSamples, truncSFS, ph,
                                  (const Polya *) polya, errTol);

	// Correct spectrum
	double sum=0.0;
    unsigned foldDim = specdim(nSamples, true);
    double truFolded[foldDim];
    for(i=0; i < foldDim; ++i) {
        j = i+1;
		unsigned two_j = 2*j;
        if(two_j < nSamples)
            truFolded[i] = 1.0/j + 1.0/(nSamples-j);
        else
            truFolded[i] = 1.0/j;
        if(i >= truncSFS)
            sum += truFolded[i];
    }

    double err = TFESpectrum_diff(tfes, foldDim, truFolded);
    assert(fabs(err) < 0.0001);

    for(i=0; i<truncSFS; ++i)
        truFolded[i] = 0.0;
    for(i=truncSFS; i<foldDim; ++i)
        truFolded[i] /= sum;

	if(verbose) {
		printf("nSamples=%u\n", nSamples);
		printf("%5s %10s %10s\n", "", "calc", "truth");
	}
	for(i=truncSFS; i < foldDim ; ++i) {
        j = i+1;
		double absErrFolded = fabs(TFESpectrum_atNdx(tfes, j) - truFolded[i]);
		int fail=0;
		if(absErrFolded > errTol) {
			fail = 1;
			ok = 0;
		}

		if(verbose){
			printf("%5u %10.5lf %10.5lf",
				   j,
				   TFESpectrum_atNdx(tfes, j),
                   truFolded[i]);
			if(fail)
				printf(" FAIL: errFolded=%lg\n", absErrFolded);
			else
				putchar('\n');
		}
	}

    TFESpectrum_free(tfes);
    Polya_free(polya);
	if(ok)
		unitTstResult("TFESpectrum", "OK");
	else
		unitTstResult("TFESpectrum", "FAIL");
}
