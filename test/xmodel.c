/**
 * @file xmodel.c
 * @author Alan R. Rogers
 * @brief Test model.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "model.h"
#include "misc.h"
#include "tokenizer.h"
#include "pophist.h"
#include "hill.h"
#include "strobeck.h"
#include "hayes.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_odeiv2.h>

int main(int argc, char **argv) {

    Model      *model;
    EpochLink  *linkedList = NULL;
    PopHist    *ph;
    int         i, twoNsmp = 20;
    double      twoN = 1e4;
    double      c = 1e-8;
    double      kb;             /* kilobases separating the two sites */
    double      u = 1e-6;
    double      s;
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;
    int         verbose = 0;
    int         ok = 1;

#ifdef NDEBUG
    eprintf("ERR@%s:%d:"
            "Unit tests must be compiled without -DNDEBUG flag\n",
            __FILE__, __LINE__);
#endif

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xmodel [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xhill [-v]\n");
    }

    model = Model_alloc("Strobeck", twoNsmp);
    ODE        *ode = ODE_new(model, odeAbsTol, odeRelTol);

    checkmem(ode, __FILE__, __LINE__);

    assert(Model_stateDim(model) == ODE_stateDim(ode));

    if(verbose) {
        printf("Initial state of model %s: ", Model_lbl(model));
        ODE_printState(ode, stdout);
        putchar('\n');
    }

    /*
     * Set up pophist in which recent epoch has t=300, N=1e2, and
     * earliest epoch has t=Inf, N=1e5,.
     */
    linkedList = EpochLink_new(linkedList, 300.0, 1e2);
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 1e5);
    ph = PopHist_fromEpochLink(linkedList);

    EpochLink_free(linkedList);
    linkedList = NULL;

    if(verbose)
        printf("%8s %8s %8s\n", "c", "LD", "LD_eq0");
    for(kb = 1.0; kb <= 300; kb += 5.0) {
        double      sep_bases = kb * 1000;

        s = ODE_ld(ode, c * sep_bases, u, ph);
        double      s0 = ODE_ldEq(ode, c * sep_bases, u, ph, 0);

        if(verbose) {
            printf("%8.6lf %8.6lf %8.6lf\n", c * sep_bases, s, s0);
        }
    }

    if(verbose) {
        printf("Final state of model via ODE_printState: ");
        ODE_printState(ode, stdout);
        putchar('\n');
        printf("And via ODE_stateVal: [");
        for(i = 0; i < ODE_stateDim(ode); ++i)
            printf(" %lg", ODE_stateVal(ode, i));
        printf("]\n");
        printf("State labels: [");
        for(i = 0; i < Model_stateDim(model); ++i)
            printf(" %s", Model_stateLbl(model, i));
        printf("]\n");
    }

    int         nbins = 10;
    double      ldvec[nbins], ldvec0[nbins], cc[nbins];
    double      hiCm = 0.3;
    double      binwidth = hiCm / nbins;

    for(i = 0; i < nbins; ++i)
        cc[i] = 0.01 * (i + 0.5) * binwidth;
    ODE_ldVec(ode, ldvec, nbins, cc, u, ph);
    ODE_ldVecEq(ode, ldvec0, nbins, cc, u, ph, 0);
    if(verbose)
        printf("%8s %8s %8s (using ODE_ldVec and ODE_ldVecEq)\n",
               "c", "LD", "LD_eq0");
    for(i = 0; verbose && i < nbins; ++i)
        printf("%8.6lg %8.6lg %8.6lg\n", cc[i], ldvec[i], ldvec0[i]);

    assert(model->stateVal(ODE_state(ode), 0) == ODE_stateVal(ode, 0));

    void       *state = Model_newState(model);
    unsigned    ydim = Model_stateDim(model);

    Strobeck_geteq(state, twoN, c, u);
    int         status = ODE_evolve(ode, state, ydim, c, u, ph,
                                    Strobeck_dydt, verbose);

    if(status) {
        printf("ODE_evolve failed\n");
        ok = 0;
    }

    Model_freeState(model, state);
    ODE_free(ode);
    Model_free(model);
    ode = NULL;

    ModelList  *ml = ModelList_alloc("Hill,Strobeck", twoNsmp);

    assert(ModelList_size(ml) == 2);
    assert(ModelList_addModel(ml, "Hayes", twoNsmp) == 1);
    assert(ModelList_size(ml) == 3);
    model = ModelList_model(ml, 1);
    assert(strcmp(Model_lbl(model), "Strobeck") == 0);
    ModelList_free(ml);
    ml = NULL;

    PopHist_free(ph);
    ph = NULL;

    if(ok)
        unitTstResult("Model", "OK");
    else
        unitTstResult("Model", "FAIL");

    return 0;
}
