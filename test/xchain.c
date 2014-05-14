/**
 * @file xchain.c
 * @author Alan R. Rogers
 * @brief Test chain.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "chain.h"
#include "pophist.h"
#include "model.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

/**
 * Log of objective function: negative of sum of differences between
 * observed and expected sigma_d^2 vectors.
 *
 * This version is a dummy with a known answer. As one moves away from
 * the global maximum at x=10, there local maxima of declining
 * magnitude at intervals of 2*Pi. Thus, the chain may get stuck at
 * x=3.7 or x=16.3. If things go well, it should find the global
 * optimum at x=10.
 *
 * @param[in] ph Current population history
 * @param[in] u mutation rate per site per generation
 * @param[in] nbins Number of values in vectors obs and c
 * @param[in] sigdsq Vector of nbins values, the observed values of
 * sigdsq.
 * @param[in] c Vector of nbins values, the recombination rates
 * associated with the values in sigdsq.
 */
double lnObjFun(PopHist * ph, double u, ODE * ode,
                int nbins, double *sigdsq, double *c) {
    double x = PopHist_twoNinv(ph, 0);
    return log((1 + cos(x - 10)) / (1 + fabs((x - 10) / 10)));
}

int main    (int argc, char **argv) {

    int verbose = 0;
    int twoNsmp = 60;
    unsigned nchains = 5;
    unsigned i;
    int rval;
    double temp = 0.8;
    unsigned rep, nreps = 100000;
    double u = 1e-8;
    double dt = 1.0;
    double dNinv = 0.2;
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned) time(NULL));

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xchain [-v]\n");
        verbose = 1; break;
    default:
        eprintf("usage: xchain [-v]\n");
    }

    double sigdsq_obs[] = {
        0.30302739, 0.20007912, 0.16397083, 0.13770996,
        0.11695362, 0.10015091, 0.08647066, 0.07487214,
        0.06469971, 0.05664014, 0.04980050, 0.04412119,
        0.03947164, 0.03532209, 0.03186997, 0.02907910,
        0.02666333, 0.02437148, 0.02251085, 0.02104841,
        0.01968679, 0.01862835, 0.01762941, 0.01669749,
        0.01601805};
    double cm[] = {
        0.00592027, 0.01799641, 0.02999884, 0.04199953,
        0.05399761, 0.06600001, 0.07799899, 0.09000018,
        0.10200012, 0.11399689, 0.12599981, 0.13799764,
        0.14999787, 0.16199999, 0.17400004, 0.18599952,
        0.19800026, 0.20999918, 0.22199915, 0.23399994,
        0.24599758, 0.25800046, 0.27000090, 0.28199899,
        0.29399900};
    unsigned nbins = sizeof(cm) / sizeof(cm[0]);
    double c[nbins];

    for(i = 0; i < nbins; ++i)
        c[i] = cm[i] * 0.01;

    EpochLink * linkedList = NULL;
    linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 3.0);
    PopHist * ph = PopHist_fromEpochLink(linkedList);
    EpochLink_free(linkedList);
    linkedList = NULL;
    Model * model = Model_alloc("Hill", twoNsmp);
    unsigned swapInterval = 1;
    pthread_t threads[nchains];
    Chain * ch, *chain;
    double odeAbsTol = 1e-7;
    double odeRelTol = 1e-3;
    chain = Chain_new(0, nchains, nreps, u, dt, dNinv, temp,
                      nbins, sigdsq_obs, c, ph, model,
                      odeAbsTol, odeRelTol, swapInterval);
    if(verbose) {
        fflush(stdout);
        fprintf(stderr, "Launching %d threads...\n", nchains);
    }
    i = 0;
    for(ch = chain; ch != NULL; ch = Chain_next(ch)) {
        rval = pthread_create(&threads[i], NULL, runChain, ch);
        if(rval)
            eprintf("ERR@%s:%d: pthread_create returned %d\n",
                    __FILE__, __LINE__, rval); ++i;
    }

    if(verbose) {
        fflush(stdout);
        fprintf(stderr, "%s: Beginning %u reps...\n", __func__, nreps);
        Chain_printHdr(chain, stdout);
    }
    for(rep = 0; rep < nreps; ++rep) {

        if(verbose)
            PRSTAT("waiting for data");
        Chain_waitForData(chain);
        if(verbose) {
            PRSTAT("got data");
            Chain_printState(chain, stdout);
        }
        Chain_setDataNeeded(chain);
        Chain_unlock(chain);
        Chain_signal(chain);
    }

    if(verbose) {
        fflush(stdout);
        fprintf(stderr, "Waiting for threads to finish...\n");
    }
    /* wait for threads to finish */
    for(i = 0; i < nchains; ++i) {
        void *status;
        rval = pthread_join(threads[i], &status);
        if(rval)
            eprintf("ERR@%s:%d: pthread_join returned %d\n",
                    __FILE__, __LINE__, rval);
        if(verbose)
            fprintf(stderr, " %2d threads have finished\n", i + 1);
    }

    if(verbose)
        fprintf(stderr, "Back from threads\n");
    PopHist * bestPh = PopHist_newEmpty(PopHist_nepoch(ph));
    double bestLnObj;
    Chain_bestFit(chain, &bestLnObj, bestPh);
    if(verbose) {
        printf("Optimum: lnObj=%g\n", bestLnObj);
        PopHist_print_comment(bestPh, "#    ", stdout);
    }
    for(ch = chain; ch != NULL; ch = Chain_next(ch)) {
        unsigned naccpt = Chain_naccepted(ch);
        unsigned nswap = Chain_nswapped(ch);
        if(verbose) {
            fprintf(stderr, "Chain %u: accepted %u/%u=%lf;"
                    " swapped %u/%u=%lf\n",
                    Chain_which(ch),
                    naccpt,
                    nreps,
                    naccpt / (double) nreps,
                    nswap, nreps, nswap / (double) nreps);
        }
    }

    PopHist_free(ph);
    PopHist_free(bestPh);
    gsl_rng_free(rng);
    Model_free(model);
    Chain_free(chain);

    unitTstResult("Chain", "OK");

    pthread_exit(NULL);
    return 0;
}
