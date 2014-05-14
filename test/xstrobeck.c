/**
 * @file xstrobeck.c
 * @author Alan R. Rogers
 * @brief Test strobeck.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strobeck.h"
#include "pophist.h"
#include "misc.h"
#include "model.h"
#include "string.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <math.h>

int main(int argc, char **argv) {

    double      c = 1e-8;

    /* number of kilobases separating the two sites */
    double      kb;

    double      u = 1e-6;
    double      relerr, relerr_tol = 0.004;
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;
    int         twoNsmp = 20;
    int         ok = 1;
    int         i, verbose = 0;
    const char *statlbl;
    size_t      stateDim = Strobeck_stateDim();
    double      ystrobeck[stateDim], y_ode[stateDim], y_eq0[stateDim];
    Model      *model = Model_alloc("Strobeck", twoNsmp);
    ODE        *ode = ODE_new(model, odeAbsTol, odeRelTol);

    assert(Model_stateDim(model) == stateDim);

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xstrobeck [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xhill [-v]\n");
    }

    /* Set up pophist */
    EpochLink  *linkedList = NULL;

    linkedList = EpochLink_new(linkedList, 300.0, 1e5);
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 1e2);
    PopHist    *ph = PopHist_fromEpochLink(linkedList);

    for(kb = 1.0; kb <= 300; kb += 5.0) {
        double      sep_bases = kb * 1000;

        Strobeck_evolveDiscrete(ystrobeck, ph, c * sep_bases, u);
        double      sigdsq_hill = Strobeck_get_sigdsq(ystrobeck, twoNsmp);

        double      sigdsq_eq0 = ODE_ldEq(ode, c * sep_bases, u, ph, 0);

        for(i = 0; i < Model_stateDim(model); ++i)
            y_eq0[i] = ODE_stateVal(ode, i);

        double      sigdsq_ode = ODE_ld(ode, c * sep_bases, u, ph);

        for(i = 0; i < Model_stateDim(model); ++i)
            y_ode[i] = ODE_stateVal(ode, i);

        if(verbose) {
            printf("%8s %8s %8s %8s %8s %8s\n",
                   "param", "c", "hill", "ode", "relerr", "eqEpoch0");
        }
        for(i = 0; i < Model_stateDim(model); ++i) {
            char        buff[10];

            snprintf(buff, sizeof(buff), "y[%d]", i);
            relerr = fabs(y_ode[i] - ystrobeck[i]) / fabs(ystrobeck[i]);
            if(verbose)
                printf("%8s %8.6f %8.6f %8.6f %8.6f %8.6f\n",
                       buff, c * sep_bases, ystrobeck[i], y_ode[i], relerr,
                       y_eq0[i]);
        }

        relerr = fabs(sigdsq_hill - sigdsq_ode) / fabs(sigdsq_hill);
        if(relerr > relerr_tol) {
            ok = 0;
            statlbl = "FAIL";
        } else
            statlbl = "OK";
        if(verbose) {
            printf("%8s %8.6f %8.6f %8.6f %8.6f %8.6f %6s\n",
                   "sigdsq", c * sep_bases, sigdsq_hill, sigdsq_ode, relerr,
                   sigdsq_eq0, statlbl);
        }
    }
    if(verbose) {
        printf("c=%lg u=%lg\n", c, u);
        PopHist_print(ph, stdout);
    }

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
    ODE_free(ode);
    EpochLink_free(linkedList);
    PopHist_free(ph);

    if(ok)
        unitTstResult("Strobeck", "OK");
    else
        unitTstResult("Strobeck", "FAIL");

    return 0;
}
