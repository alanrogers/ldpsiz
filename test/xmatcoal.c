/**
   @file xmatcoal.c
   @brief Test matcoal.c
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "matcoal.h"
#include "pophist.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    unsigned i, j, nSamples=30;
    int verbose = 0, ok=1;
    double err, maxerr, errTol = 1e-6;

    for(j=1; j<argc; ++j) {
        if(strcmp(argv[j], "-v") == 0)
            verbose = 1;
        else if(strcmp(argv[j], "-n") == 0) {
            if(++j >= argc) {
                fprintf(stderr,"missing arg for -n\n");
                exit(1);
            }

            errno = 0;
            nSamples = strtoul(argv[j], 0, 10);
            if(errno) {
                perror("strtoul");
                exit(EXIT_FAILURE);
            }
        }
    }
	
    EpochLink  *linkedList = EpochLink_new(NULL, HUGE_VAL, 1.0);
    PopHist *ph = PopHist_fromEpochLink(linkedList);
    checkmem(ph, __FILE__, __LINE__);
    PopHist_sanityCheck(ph, __FILE__, __LINE__);

    if(verbose)
        PopHist_print(ph, stdout);

    long double x[nSamples];
    long double betavec[nSamples];
    for(i=0; i<nSamples; ++i)
        betavec[i] = MatCoal_beta(i);

    /* construct initial prob vector */
    memset(x, 0, (nSamples-1)*sizeof(x[0]));
    x[nSamples-1] = 1.0;  /* initially all prob is at end */

    /* project 0 units into future */
    MatCoal_project(nSamples, x, 0.0, betavec, errTol);

    /* check result */
    if(verbose)
        printf("t=0:\n");
    maxerr = 0.0;
    for(i=0; i<nSamples-1; ++i) {
        err = fabs(x[i] - 0.0);
        if(err > maxerr)
            maxerr = err;
        if(verbose) {
            fprintf(stderr, "x[%d]=%Lf; should=%Lf; err=%lf\n",
                    i, x[i], 0.0L, err);
        }
    }
    err = fabs(x[nSamples-1] - 1.0);
    if(err > maxerr)
        maxerr = err;
    if(verbose) {
        fprintf(stderr, "x[%u]=%Lf; should=%Lf; err=%lf\n",
                nSamples-1, x[nSamples-1], 1.0L, err);
    }
    if(maxerr > errTol)
        ok = 0;
		
    /* re-construct initial prob vector */
    memset(x, 0, (nSamples-1)*sizeof(x[0]));
    x[nSamples-1] = 1.0;  /* initially all prob is at end */

    /* project 1 units into future */
    MatCoal_project(nSamples, x, 1.0, betavec, errTol);

    if(verbose) {
        /* print projection */
        printf("t=1:\n");
        for(i=0; i<nSamples; ++i)
            printf("  x[%d] = %Lf\n", i+1, x[i]);
    }

    /* project from t=1 to t=2 */
    MatCoal_project(nSamples, x, 1.0, betavec, errTol);

    long double y[nSamples];

    /* re-construct initial prob vector */
    memset(y, 0, (nSamples-1)*sizeof(y[0]));
    y[nSamples-1] = 1.0L;  /* initially all prob is at end */

    /* project from 0 to 2 in 1 step */
    MatCoal_project(nSamples, y, 2.0, betavec, errTol);

    /* x and y should be the same (both are x(2)) */
    maxerr = 0.0;
    if(verbose)
        printf(" t=2: %15s %15s\n", "2-step", "1-step");
    for(i=0; i < nSamples; ++i) {
        err = fabs(x[i] - y[i]);
        if(err > maxerr)
            maxerr = err;
        if(verbose)
            printf("      %15.8Lg %15.8Lg\n", x[i], y[i]);
    }
    if(verbose)
        printf("t=2: maxerr=%lf\n", maxerr);
    if(maxerr > errTol)
        ok = 0;

    /* project from t=2 to t=3 */
    MatCoal_project(nSamples, x, 1.0, betavec, errTol);

    /* re-construct initial prob vector */
    memset(y, 0, (nSamples-1)*sizeof(y[0]));
    y[nSamples-1] = 1.0;  /* initially all prob is at end */

    /* project from 0 to 3 in 1 step */
    MatCoal_project(nSamples, y, 3.0, betavec, errTol);

    /* x and y should be the same (both are x(3)) */
    maxerr = 0.0;
    if(verbose)
        printf(" t=3: %15s %15s\n", "3-step", "1-step");

    for(i=0; i < nSamples; ++i) {
        err = fabs(x[i] - y[i]);
        if(err > maxerr)
            maxerr = err;
        if(verbose)
            printf("%5d %15.8Lg %15.8Lg\n", i+1, x[i], y[i]);
    }
    if(verbose)
        printf("t=3: maxerr=%lg\n", maxerr);
    if(maxerr > errTol)
        ok = 0;

    unitTstResult("MatCoal_project", ok ? "OK" : "FAIL");

    ok = 1;

    /* Test project(Matrix, vector) */
    unsigned nTimes=5;
    double m[nTimes][nSamples];
    
    /* row 0 is (0,...,0,1) */
    for(i=0; i < nSamples-1; ++i)
        m[0][i] = 0.0;
    m[0][nSamples-1] = 1.0;

    double tvec[nTimes];

    for(i=0; i < nTimes; ++i)
        tvec[i] = i;
    tvec[nTimes-1] += 10.0;

    /*
     * Project initial vector tvec forward to a succession of times.
     *
     * tvec: an array of nTimes doubles, giving the time values.
     *
     * m: a nSamplesXnTimes matrix of doubles. On return, m[i][j] contains
     * the probability that the coalescent contains j+1 distinct
     * lineages at time tvec[i].
     *
     */
    MatCoal_project_multi(nTimes, nSamples, m, tvec, ph, errTol);

    if(verbose) {

        printf("Cols show probabilities at times given in time vector:\n");
        char buff[20];
        printf("%4s", "n");
        for(j=0; j < nTimes; ++j) {
            snprintf(buff, sizeof(buff), "t=%lg", tvec[j]);
            printf(" %15s", buff);
        }
        putchar('\n');
        for(i=0; i<nSamples; ++i) {
            printf("%4u", i+1);
            for(j=0; j < nTimes; ++j)
                printf(" %15.8g", m[j][i]);
            putchar('\n');
        }
    }

    unitTstResult("MatCoal_project_multi", "OK");
    ok = 1;

    /* Test integrate */
    double z[nSamples];

    MatCoal_integrate(nSamples, z, ph, errTol);

    for(i=2; i<nSamples; ++i) {
        double truth = 2.0/(i*(i-1));
        err = fabs( z[i-1] - truth );
        if(err >= errTol) {
            printf("z[%d]=%lf", i-1, z[i-1]);
            printf("  truth=%lf", truth);
            printf("  diff=%lf\n", fabs( z[i-1] - truth ));
            ok = 0;
        }
    }

    if(verbose) {
        for(i=0; i<nSamples; ++i)
            printf("E[length of interval %d] = %6lg; should be %lg\n",
                   i+1, z[i], 2.0/(i*(i+1)));
    }

    unitTstResult("MatCoal_integrate", ok ? "OK" : "FAIL");
    return 0;
}





