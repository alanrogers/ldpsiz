/**
   @file xespectrum.c
   @brief Test class ESpectrum
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "espectrum.h"
#include "pophist.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

int main(int argc, char **argv) {
    unsigned nSamples=10u;
    int verbose = 0, ok=1;
	unsigned i;
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
	
	// pophist has 1 epoch with pop size twoN=1.
    EpochLink  *linkedList = EpochLink_new(NULL, HUGE_VAL, 1.0);
    PopHist *ph = PopHist_fromEpochLink(linkedList);
    ESpectrum *es = ESpectrum_new(nSamples, ph, errTol);

	// Correct spectrum
	double truUnfolded[nSamples];
	double sum=0.0;
	for(i=1; i < nSamples; ++i)
		sum += truUnfolded[i] = 1.0/i;
	for(i=1; i < nSamples; ++i)
		truUnfolded[i] /= sum;

	if(verbose) {
		printf("nSamples=%u\n", nSamples);
        printf("%5s %21s %21s\n", "", "Unfolded", "Folded");
		printf("%5s %10s %10s %10s %10s\n", "",
			   "calc", "truth", "calc", "truth");
	}
	for(i=1; i < nSamples; ++i) {
		double truFolded;
		unsigned two_i = 2*i;
		if(two_i < nSamples)
			truFolded = truUnfolded[i]+truUnfolded[nSamples-i];
		else if(two_i == nSamples)
			truFolded = truUnfolded[i];
		else
			truFolded = 0.0;
		double absErrFolded = fabs(ESpectrum_folded(es, i) - truFolded);
		double absErrUnfolded = fabs(ESpectrum_unfolded(es, i)
									 - truUnfolded[i]);
		int fail=0;
		if(absErrFolded > errTol
		   || absErrUnfolded > errTol) {
			fail = 1;
			ok = 0;
		}

		if(verbose){
			printf("%5u %10.5lf %10.5lf %10.5lf %10.5lf",
				   i,
				   ESpectrum_unfolded(es, i),
				   truUnfolded[i],
				   ESpectrum_folded(es, i),
                   truFolded);
			if(fail)
				printf(" FAIL: errUnfolded=%lg errFolded=%lg\n",
					   absErrUnfolded, absErrFolded);
			else
				putchar('\n');
		}
	}

	if(ok)
		unitTstResult("ESpectrum", "OK");
	else
		unitTstResult("ESpectrum", "FAIL");
}





