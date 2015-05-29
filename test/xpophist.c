/**
 * @file xpophist.c
 * @author Alan R. Rogers
 * @brief Test pophist.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "pophist.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>

int main(int argc, char **argv) {

    int         i, verbose = 0;
    double      t, t2;
    double      twoN, twoN1, twoN2, twoNinv1, twoNinv2;
    char        buff[100];
    PopHist    *ph = NULL;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xpophist [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xpophist [-v]\n");
    }

    EpochLink_test();
    unitTstResult("EpochLink", "OK");

    int         nepoch = 2;

    ph = PopHist_newEmpty(nepoch);
    assert(nepoch == PopHist_nepoch(ph));
    if(verbose) {
        PopHist_print(ph, stdout);
        PopHist_print_comment(ph, "comment>", stdout);
    }
    for(i = 0; i < nepoch; ++i) {
        assert(0.0 == PopHist_twoNinv(ph, i));
        t = PopHist_duration(ph, i);
        if(i == nepoch - 1) {
            assert(isinf(t) && t > 0.0);
        } else
            assert(0.0 == t);
    }
    if(verbose)
        unitTstResult("PopHist_newEmpty", "OK");

    PopHist_free(ph);
    ph = NULL;

    EpochLink *el = NULL;
    el = EpochLink_new(el, 50.0, 1000.0);
    el = EpochLink_new(el, 50.0, 1000.0);
    el = EpochLink_new(el, 10.0, 1000.0);   /* duration 10 becomes inf */
    ph = PopHist_fromEpochLink(el);
    nepoch = 3;
    EpochLink_free(el);
    el = NULL;

    assert(3 == PopHist_nepoch(ph));
    twoN = PopHist_twoN(ph, 0);
    assert(1000.0 == twoN);
    twoNinv1 = PopHist_twoNinv(ph, 0);
    assert(1.0 / 1000.0 == twoNinv1);
    assert(50.0 == PopHist_duration(ph, 0));
    twoN = PopHist_twoN(ph, 1);
    assert(1000.0 == twoN);
    twoNinv1 = PopHist_twoNinv(ph, 1);
    assert(1.0 / 1000.0 == twoNinv1);
    assert(50.0 == PopHist_duration(ph, 1));
    twoN = PopHist_twoN(ph, 2);
    assert(1000.0 == twoN);
    assert(1.0 / 1000.0 == PopHist_twoNinv(ph, 2));
    t = PopHist_duration(ph, 2);
    assert(isinf(t) && t > 0.0);

    if(verbose)
        unitTstResult("PopHist_fromEpochLink", "OK");

    PopHist_setDuration(ph, 0, 222.0);
    assert(222.0 == PopHist_duration(ph, 0));

    if(verbose)
        unitTstResult("PopHist_setDuration", "OK");

    int         npar = PopHist_nParams(ph);

    assert(2 * nepoch - 1 == npar);
    if(verbose)
        unitTstResult("PopHist_nParams", "OK");
    if(verbose) {
        for(i = 0; i < npar; ++i) {
            PopHist_paramName(ph, buff, sizeof(buff), i);
            printf("parameter %d: %s = %lg\n",
                   i, buff, PopHist_paramValue(ph, i));
        }
        unitTstResult("PopHist_paramName", "OK");
        unitTstResult("PopHist_paramValue", "OK");
    }

    gsl_vector *state = gsl_vector_alloc((unsigned) npar);

    PopHist    *ph2 = PopHist_newEmpty(nepoch);
    PopHist_to_vector(state, ph);
    vector_to_PopHist(ph2, state);

    PopHist    *ph3 = PopHist_newEmpty(nepoch);
    double      stateArray[npar];
    PopHist_to_C_array(npar, stateArray, ph);
    C_array_to_PopHist(ph3, npar, stateArray);

    assert(0.0 == PopHist_distance(ph, ph));
    assert(0.0 == PopHist_distance(ph, ph2));
    assert(0.0 == PopHist_distance(ph, ph3));

    for(i = 0; i < nepoch; ++i) {
        twoNinv1 = PopHist_twoNinv(ph, i);
        twoNinv2 = PopHist_twoNinv(ph2, i);
        assert(twoNinv1 == twoNinv2);
        twoNinv2 = PopHist_twoNinv(ph3, i);
        assert(twoNinv1 == twoNinv2);
        if(i < nepoch - 1) {
            assert(PopHist_duration(ph, i) == PopHist_duration(ph2, i));
            assert(PopHist_duration(ph, i) == PopHist_duration(ph3, i));
        } else {
            t = PopHist_duration(ph, i);
            assert(isinf(t) && t > 0.0);
            t = PopHist_duration(ph2, i);
            assert(isinf(t) && t > 0.0);
            t = PopHist_duration(ph3, i);
            assert(isinf(t) && t > 0.0);
        }
    }
    if(verbose) {
        unitTstResult("PopHist_to_vector", "OK");
        unitTstResult("PopHist_to_C_array", "OK");
        unitTstResult("vector_to_PopHist", "OK");
        unitTstResult("C_array_to_PopHist", "OK");
    }

    PopHist_setDuration(ph, 0, 321.0);
    PopHist_setTwoN(ph, 0, 432.0);
    PopHist_setTwoN(ph, 1, 543.0);

    assert(0.0 < PopHist_distance(ph, ph2));

    PopHist_copy(ph2, ph);
    assert(0.0 == PopHist_distance(ph, ph2));
    for(i = 0; i < nepoch; ++i) {
        twoN1 = PopHist_twoN(ph, i);
        twoN2 = PopHist_twoN(ph2, i);
        assert(twoN1 == twoN2);
        t = PopHist_duration(ph, i);
        t2 = PopHist_duration(ph2, i);
        if(i < nepoch - 1) {
            assert(t == t2);
        } else {
            assert(isinf(t) && isinf(t2) && t > 0.0 && t2 > 0.0);
        }
    }
    if(verbose)
        unitTstResult("PopHist_copy", "OK");

    PopHist_free(ph2);
    ph2 = NULL;

    PopHist_setDuration(ph, 0, 888.0);
    PopHist_setTwoN(ph, 0, 777.0);
    PopHist_setTwoN(ph, 1, 666.0);

    ph2 = PopHist_dup(ph);
    assert(0.0 == PopHist_distance(ph, ph2));
    for(i = 0; i < nepoch; ++i) {
        twoN1 = PopHist_twoN(ph, i);
        twoN2 = PopHist_twoN(ph2, i);
        assert(twoN1 == twoN2);
        if(i < nepoch - 1) {
            assert(PopHist_duration(ph, i) == PopHist_duration(ph2, i));
        } else {
            t = PopHist_duration(ph, i);
            t2 = PopHist_duration(ph2, i);
            assert(isinf(t) && isinf(t2) && t > 0.0 && t2 > 0.0);
        }
    }
    if(verbose)
        unitTstResult("PopHist_dup", "OK");

    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

    gsl_rng_set(rng, (unsigned) time(NULL));

    double      dt = 1.0, dNinv = 0.03;

    PopHist_perturb(ph2, ph, dt, dNinv, rng);
    assert(0.0 < PopHist_distance(ph, ph2));
    for(i = 0; i < nepoch; ++i) {
        assert(PopHist_twoN(ph, i) != PopHist_twoN(ph2, i));
#ifndef PERTURB_GAUSSIAN
#ifndef PERTURB_TDIST
        assert(2 * dNinv >= fabs(PopHist_twoNinv(ph, i)
                                 - PopHist_twoNinv(ph2, i)));
#endif
#endif
        if(i < nepoch - 1) {
            assert(PopHist_duration(ph, i) != PopHist_duration(ph2, i));
#ifndef PERTURB_GAUSSIAN
#ifndef PERTURB_TDIST
            assert(2 * dt >= fabs(PopHist_duration(ph, i)
                                  - PopHist_duration(ph2, i)));
#endif
#endif
        } else {
            t = PopHist_duration(ph, i);
            t2 = PopHist_duration(ph2, i);
            assert(isinf(t) && isinf(t2) && t > 0.0 && t2 > 0.0);
        }
    }
    if(verbose)
        unitTstResult("PopHist_perturb", "OK");

    if(verbose)
        unitTstResult("PopHist_distance", "OK");

    gsl_rng_free(rng);
    gsl_vector_free(state);
    PopHist_free(ph);
    PopHist_free(ph2);

    unitTstResult("PopHist", "OK");

    return 0;
}
