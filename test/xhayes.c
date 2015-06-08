/**
 * @file xhayes.c
 * @author Alan R. Rogers
 * @brief Test hayes.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "hayes.h"
#include "model.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <string.h>

int main(int argc, char **argv) {

    EpochLink  *linkedList = NULL;
    PopHist    *ph;
    double      c = 1e-8;
    double      u = 1e-6;
    int         twoNsmp = 0;
    int         ok = 1;
    int         verbose = 0;
    double      kb = 200.0;     /* kilobases */

    Model      *model = Model_alloc("Hayes", twoNsmp);

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xhill [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xhill [-v]\n");
    }

    /* Set up pophist */
    linkedList = EpochLink_new(linkedList, 300.0, 1e5);
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 1e2);
    ph = PopHist_fromEpochLink(linkedList);

    if(verbose) {
#ifdef HAYES_MUTATION_ADJUSTMENT
        printf("Adjusting for mutation: E[rsq] = 1/(4Nc+2)\n");
#else
        printf("Not adjusting for mutation: E[rsq] = 1/(4Nc+1)\n");
#endif
        printf("%8s %8s %9s\n", "c", "Sved", "Sved");
    }
    for(kb = 1.0; kb <= 300; kb += 20.0) {
        double      sep_bases = kb * 1000;

        if(verbose) {
            printf("%8.6f %8.6f %9.6f\n",
                   c * sep_bases,
                   Hayes_rsq(NULL, c * sep_bases, u, ph, 0),
                   Hayes_rsqNoODE( c * sep_bases, u, ph, twoNsmp, NULL)
                );
        }
    }
    if(verbose) {
        printf("c=%lg u=%lg twoNsmp=%d\n", c, u, twoNsmp);
        PopHist_print(ph, stdout);
    }

    assert(Model_stateDim(model) == 0);

#if 0
    /*
     * Neither of these should ever be called. If they are called by
     * mistake, they abort. Turn this code on if you want to make
     * sure they abort as they should.
     */
    Model_stateLbl(model, 0);   /* aborts */
    Model_state(model, 0);      /* aborts */
#endif

    int         i, nbins = 10;
    double      rsq[nbins], cc[nbins], eq0[nbins], eq1[nbins];
    double      hiCm = 0.3;
    double      width = hiCm / nbins;

    for(i = 0; i < nbins; ++i) {
        cc[i] = 0.01 * (i + 0.5) * width;
        rsq[i] = Hayes_rsq(NULL, cc[i], u, ph, twoNsmp);
        eq0[i] = Hayes_rsqEq(cc[i], u, ph, 0, twoNsmp, NULL);
        eq1[i] = Hayes_rsqEq(cc[i], u, ph, 1, twoNsmp, NULL);
    }
    if(verbose) {
        printf("Result from Hayes_rsq and Hayes_rsqEq:\n");
        printf("%8s %10s %10s %10s\n", "c", "rsq", "eq0", "eq1");
        for(i = 0; i < nbins; ++i)
            printf("%8.6lg %10.6lg %10.6lg %10.6lg\n",
                   cc[i], rsq[i], eq0[i], eq1[i]);
    }

    Model_free(model);

    if(ok)
        printf("Hayes OK\n");
    else {
        printf("Hayes FAIL\n");
        exit(EXIT_FAILURE);
    }

    return 0;
}
