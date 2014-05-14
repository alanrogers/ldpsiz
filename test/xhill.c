/**
 * @file xhill.c
 * @brief Test hill.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "hill.h"
#include "misc.h"
#include "pophist.h"
#include "model.h"
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <float.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char **argv) {

    EpochLink  *linkedList = NULL;
    PopHist    *ph;
    int         nbins = 20;
    size_t      stateDim = Hill_stateDim();

    double      c, cm;
    double      hi_cm = 0.3;
    double      lo_cm = hi_cm / (2 * nbins);
    double      step_cm = (hi_cm - lo_cm) / (nbins - 1);

    double      u = 1e-6;
    double      relerr, relerr_tol = 0.004;
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;
    int         twoNsmp = 20;
    int         ok = 1;
    int         i, verbose = 0;
    const char *statlbl;
    double      y_hill[stateDim], y_ode[stateDim], y_eq0[stateDim];
    Model      *model = Model_alloc("Hill", twoNsmp);
    ODE        *ode = ODE_new(model, odeAbsTol, odeRelTol);

    assert(Model_stateDim(model) == stateDim);

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

#if 0
    linkedList = EpochLink_new(linkedList, 300.0, 1e5);
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 1e2);
#else
    linkedList = EpochLink_new(linkedList, 731.076, 322.111);
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 845.946);
#endif
    ph = PopHist_fromEpochLink(linkedList);

    /*    for(kb=1.0; kb <= 300; kb += 5.0) { */
    for(cm = lo_cm; cm <= hi_cm; cm += step_cm) {
        c = cm * 0.01;

        Hill_evolveDiscrete(y_hill, ph, c, u);
        double      sigdsq_hill = Hill_get_sigdsq(y_hill, twoNsmp);

        double      sigdsq_eq0 = ODE_ldEq(ode, c, u, ph, 0);

        for(i = 0; i < Model_stateDim(model); ++i)
            y_eq0[i] = ODE_stateVal(ode, i);

        double      sigdsq_ode = ODE_ld(ode, c, u, ph);

        for(i = 0; i < Model_stateDim(model); ++i)
            y_ode[i] = ODE_stateVal(ode, i);

        if(verbose) {
            printf("%9s %9s %9s %9s %9s %9s\n",
                   "param", "cm", "hill", "ode", "relerr", "eqEpoch0");
        }
        for(i = 0; i < Model_stateDim(model); ++i) {
            char        buff[10];

            snprintf(buff, sizeof(buff), "y[%d]", i);
            relerr = fabs(y_ode[i] - y_hill[i]) / fabs(y_hill[i]);
            if(verbose)
                printf("%9s %9.6f %9.6f %9.6f %9.6f %9.6f\n",
                       Model_stateLbl(model, i),
                       /*buff, */
                       cm, y_hill[i], y_ode[i], relerr, y_eq0[i]);
        }

        relerr = fabs(sigdsq_hill - sigdsq_ode) / fabs(sigdsq_hill);
        if(relerr > relerr_tol) {
            ok = 0;
            statlbl = "FAIL";
        } else
            statlbl = "OK";
        if(verbose) {
            printf("%9s %9.6f %9.6f %9.6f %9.6f %9.6f %6s\n",
                   "sigdsq", cm, sigdsq_hill, sigdsq_ode, relerr,
                   sigdsq_eq0, statlbl);
        }
    }
    if(verbose)
        PopHist_print(ph, stdout);

    if(verbose)
        printf("State labels:");
    for(i = 0; i < Model_stateDim(model); ++i) {
        const char *lbl = Model_stateLbl(model, i);

        if(verbose)
            printf(" %s", lbl);
    }
    if(verbose)
        putchar('\n');

    Model_free(model);

    double      twoN;

    /* test getDRM */
    for(twoN = 10.0; twoN < 1e8; twoN *= 100.0) {
        for(c = 0.1; c >= 1e-8; c *= -0.01) {
            for(u = 0.01; u >= 1e-9; u *= -0.01) {
                if(0 != Hill_test_getDRM(twoN, c, u)) {
                    printf("%s:%d:%s: Hill_test_getDRM(twoN=%lg, c=%lg, u=%lg) FAIL\n",
                           __FILE__,__LINE__,__func__, twoN, c, u);
                    ok = 0;
                }
            }
        }
    }

    time_t      currtime = time(NULL);
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    double      x1[stateDim];

    gsl_rng_set(rng, (unsigned) currtime);

    /* test onestep */
    for(twoN = 10.0; twoN < 1e8; twoN *= 100.0) {
        for(c = 0.1; c >= 1e-8; c *= -0.01) {
            for(u = 0.01; u >= 1e-9; u *= -0.01) {
                for(i = 0; i < stateDim; ++i)
                    x1[i] = gsl_ran_flat(rng, 0.0, 1.0);
                if(0 != Hill_test_onestep(x1, stateDim, twoN, c, u)) {
                    printf("%s:%d:%s: Hill_test_onestep(twoN=%lg, c=%lg, u=%lg) FAIL\n",
                           __FILE__,__LINE__,__func__, twoN, c, u);
                    ok = 0;
                }
            }
        }
    }

    ODE_free(ode);
    gsl_rng_free(rng);
    EpochLink_free(linkedList);
    PopHist_free(ph);

    if(ok)
        unitTstResult("Hill", "OK");
    else {
        unitTstResult("Hill", "FAIL");
        exit(EXIT_FAILURE);
    }

    return 0;
}
