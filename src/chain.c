/**
 * @file chain.c
 * @author Alan R. Rogers
 * @brief A linked list of Markov chains, used in Markov-coupled Markov
 * chain Monte Carlo (MCMCMC). 
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include "pophist.h"
#include "model.h"
#include "chain.h"
#include "misc.h"

pthread_mutex_t stdout_mutex = PTHREAD_MUTEX_INITIALIZER;

/** The state of a single Markov chain. */
typedef struct ChainState {
    /*
     * ph[0] and ph[1] are used internally by runChain. ph[2] is used for
     * communication with other chains. Only ph[2] must be locked before
     * use.
     *
     * lnObj[0] is used internally. lnObj[1] is for communication btw
     * chains.
     */
    PopHist    *ph[3];          /**< array of 3 pointers to PopHist objects */
    unsigned    accept;         /**< 1 if state changed; 0 otherwise */
    double      lnObj[2];       /**< measure of fit to data */
} ChainState;

/** A linked list of Markov chains, which differ in the extent to
 *  which the objective function has been artificially flattened.
 */
struct Chain {
    struct Chain *next; /**< points to next Chain in linked list */

    /* invariant parameters */
    unsigned    which;     /**< index of chain in linked list */
    double      flatness;  /**< how much this chain is flattened */
    double      u;         /**< mutation rate */
    double      dt;        /**< size of perturbations in epoch lengths*/
    double      dNinv;     /**< size of perturbations in 1/2N */
    ODE        *ode;
    unsigned    nreps;     /**< number of iterations of chain */
    unsigned    swapInterval;   /**< how often to try swapping */

    /* observed data. Arrays are not locally owned */
    unsigned    nbins;   /**< dimension of arrays */
    double     *sigdsq;  /**< array of \f$\sigma_d^2$\f values*/
    double     *c;       /**< array of recombination rates */

    /* Communication between threads */
    unsigned    dataAvailable;  /**< 0:unavailable; 1:available */
    pthread_mutex_t lock;       /**< lock chain before changing state */
    pthread_cond_t signal;      /**< interthread communication */

    /*
     * Monitor behavior of chain. These variables are changed by one
     * thread only and are not consulted until all threads are
     * finished.  No need to lock.
     */
    unsigned    naccpt, nswap;

    /*
     * The following state variables change with each iteration and
     * must be locked before access.
     */
    ChainState *state;

    double      bestLnObj;
    PopHist    *bestPh;
};

ChainState *ChainState_new(PopHist * ph);
void        ChainState_free(ChainState * cs);
void        Chain_publish(Chain * chain);

/**
 * Allocate and initialize a new ChainState structure. Occupies a
 * single block of memory, using the "struct hack" of C programming.
 */
ChainState *ChainState_new(PopHist * ph) {
    int         nepoch = PopHist_nepoch(ph);
    size_t      size = sizeof(ChainState) + 3 * PopHist_calcSize(nepoch);
    ChainState *cs = malloc(size);

    memset(cs, 0, size);

    /*
     * ph[0] and ph[1] reside just after ChainState structure on same block
     * of contiguous memory.
     */
    size_t      phsize = PopHist_calcSize(nepoch);

    cs->ph[0] = (PopHist *) & cs[1];
    cs->ph[1] = (PopHist *) ((size_t) cs->ph[0] + phsize);
    cs->ph[2] = (PopHist *) ((size_t) cs->ph[0] + 2 * phsize);

    assert(cs->ph[0]);
    assert(cs->ph[1]);
    assert(cs->ph[2]);

    /* initialize PopHist structures */
    PopHist_init(cs->ph[0], nepoch, phsize);
    PopHist_init(cs->ph[1], nepoch, phsize);
    PopHist_init(cs->ph[2], nepoch, phsize);
    PopHist_copy(cs->ph[0], ph);

    return cs;
}

/** De-allocate a ChainState object */
void ChainState_free(ChainState * cs) {
    free(cs);
}

/** Make data available to other threads */
void Chain_publish(Chain * chain) {
    PopHist_copy(chain->state->ph[2], chain->state->ph[0]);
    chain->state->lnObj[1] = chain->state->lnObj[0];
    Chain_setDataAvailable(chain);
}

/** Return pointer to next Chain in linked list */
Chain      *Chain_next(Chain * chain) {
    return chain->next;
}

/** Number of proposals accepted so far by this chain */
unsigned Chain_naccepted(Chain * chain) {
    return chain->naccpt;
}

/** Number of times this chain has swapped states with the next */
unsigned Chain_nswapped(Chain * chain) {
    return chain->nswap;
}

/** Return the index of this chain */
unsigned Chain_which(Chain * chain) {
    return chain->which;
}

/** Allocate memory for a Chain */
Chain      *Chain_new(unsigned which, unsigned nChains, unsigned nreps,
                      double u, double dt, double dNinv, double temp,
                      unsigned nbins, double *sigdsq, double *c, PopHist * ph,
                      Model * model, double odeAbsTol, double odeRelTol,
                      unsigned swapInterval) {

    if(which == nChains)
        return NULL;

    Chain      *chain = malloc(sizeof(Chain));

    chain->next = Chain_new(which + 1, nChains, nreps, u, dt, dNinv, temp,
                            nbins, sigdsq, c, ph, model, odeAbsTol,
                            odeRelTol, swapInterval);

    chain->which = which;
    chain->flatness = 1.0 / (1.0 + which * temp);
    chain->u = u;
    chain->dt = dt;
    chain->dNinv = dNinv;
    chain->ode = ODE_new(model, odeAbsTol, odeRelTol);
    chain->nreps = nreps;
    chain->swapInterval = swapInterval;

    /* externally owned: not copied */
    chain->nbins = nbins;
    chain->sigdsq = sigdsq;
    chain->c = c;

    int         i;

    if((i = pthread_mutex_init(&chain->lock, NULL)))
        eprintf("ERR@%s:%d: bad return (%d) from pthread_mutex_init",
                __FILE__, __LINE__, i);

    if((i = pthread_cond_init(&chain->signal, NULL)))
        eprintf("ERR@%s:%d: couldn't init rsignal (%d)",
                __FILE__, __LINE__, i);

    chain->dataAvailable = 0;
    chain->naccpt = chain->nswap = 0;
    chain->state = ChainState_new(ph);
    chain->state->lnObj[0] = lnObjFun(chain->state->ph[0], chain->u,
                                      chain->ode, chain->nbins,
                                      chain->sigdsq, chain->c);

    chain->bestLnObj = -DBL_MAX;
    chain->bestPh = PopHist_dup(chain->state->ph[0]);

    return chain;
}

void Chain_sanityCheck(Chain * chain, const char *file, int line) {
#ifndef NDEBUG
    if(chain == NULL)
        return;

    REQUIRE(chain->which < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->flatness >= 0.0, __FILE__, __LINE__);
    REQUIRE(chain->u > 0.0, __FILE__, __LINE__);
    REQUIRE(chain->dt >= 0.0, __FILE__, __LINE__);
    REQUIRE(chain->dNinv >= 0.0, __FILE__, __LINE__);
    REQUIRE(chain->ode != NULL, __FILE__, __LINE__);
    REQUIRE(chain->swapInterval < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->nbins < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->sigdsq != NULL, __FILE__, __LINE__);
    REQUIRE(chain->c != NULL, __FILE__, __LINE__);
    REQUIRE(chain->dataAvailable < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->naccpt < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->nswap < UINT_MAX / 2u, __FILE__, __LINE__);
    REQUIRE(chain->state != NULL, __FILE__, __LINE__);
    REQUIRE(chain->bestPh != NULL, __FILE__, __LINE__);

    Chain_sanityCheck(chain->next, __FILE__, __LINE__);
#endif
}

/** Free memory allocated by linked list of chains */
void Chain_free(Chain * chain) {
    if(chain == NULL)
        return;
    Chain_free(chain->next);
    chain->next = NULL;

    pthread_cond_destroy(&chain->signal);
    pthread_mutex_destroy(&chain->lock);
#if 0
    fprintf(stderr, "%s:%d Chain_state %p freed by Chain %p\n",
            __FILE__, __LINE__, chain->state, chain);
#endif

    ChainState_free(chain->state);
    PopHist_free(chain->bestPh);
    ODE_free(chain->ode);

    free(chain);
}

/** Print machine address of state of each chain in linked list */
void Chain_prStateAddr(const Chain * chain, int ndx,
                       const char *file, int line) {
    if(chain == NULL)
        return;
    fprintf(stderr, "%s:%d: %p=state[%p] in link %d\n",
            file, line, chain->state, chain, ndx);
    Chain_prStateAddr(chain->next, ndx + 1, file, line);
}

/** Return 1 if chain's state changed; 0 otherwise */
int Chain_accepted(const Chain * c) {
    return (0 != c->state->accept);
}

/**
 * Swap states of chains c and c->next. Lock both chains before calling
 * this function.
 */
void Chain_swapState(Chain * c) {

    assert(c->next);

    ChainState *hold;

    hold = c->state;
    c->state = c->next->state;
    c->next->state = hold;
    ++c->nswap;
}

/**
 * Wait until thread has completed another iteration and made the
 * resulting data available. On return, chain will be locked.
 */
void Chain_waitForData(Chain * chain) {
    int         status;

    pthread_mutex_lock(&chain->lock);
    while(chain->dataAvailable == 0) {
        status = pthread_cond_wait(&chain->signal, &chain->lock);
        if(status)
            ERR(status, "wait for dataAvailable");
    }
    return;
}

/**
 * Wait until main thread has finished with data from last iteration.
 * On return, the chain will be locked.
 */
void Chain_waitUntilDataNeeded(Chain * chain) {
    pthread_mutex_lock(&chain->lock);
    while(chain->dataAvailable == 1)
        pthread_cond_wait(&chain->signal, &chain->lock);
}

/**
 * Indicate that data are available. Chain should be locked before
 * calling this function.
 */
void Chain_setDataAvailable(Chain * chain) {
    chain->dataAvailable = 1;
}

/**
 * Indicate that data are needed. Chain should be locked before
 * calling this function.
 */
void Chain_setDataNeeded(Chain * chain) {
    chain->dataAvailable = 0;
}

/**
 * Send signal via condition variable.
 */
void Chain_signal(Chain * chain) {
    pthread_cond_signal(&chain->signal);
}

/** Lock chain */
void Chain_lock(Chain * chain) {
    pthread_mutex_lock(&chain->lock);
    return;
}

/** Unlock chain */
void Chain_unlock(Chain * chain) {
    pthread_mutex_unlock(&chain->lock);
    return;
}

/** Lock stdout */
void Chain_lockStdout(void) {
    pthread_mutex_lock(&stdout_mutex);
    return;
}

/** Unock stdout */
void Chain_unlockStdout(void) {
    pthread_mutex_unlock(&stdout_mutex);
    return;
}

/**
 * Log of objective function. Does not calculate function. Just returns
 * current stored value.
 */
double Chain_lnObj(Chain * chain) {
    return chain->state->lnObj[1];
}

/** Return flatness parameter of chain. */
double Chain_flatness(Chain * chain) {
    return chain->flatness;
}

/** Print entire chain */
void Chain_printFull(Chain * chain, FILE * fp) {
    unsigned    i;

    fprintf(fp, "Chain %p:\n", chain);
    fprintf(fp, "  nreps=%u nbins=%u accept=%u u=%g dt=%g\n"
            "  dNinv=%g flatness=%g datAvail=%u\n",
            chain->nreps, chain->nbins, chain->state->accept,
            chain->u, chain->dt,
            chain->dNinv, chain->flatness, chain->dataAvailable);
    fprintf(fp, "  %12s %8s\n", "c", "sigdsq");
    for(i = 0; i < chain->nbins; ++i)
        fprintf(fp, "  %12.9lf %8.6lf\n", chain->c[i], chain->sigdsq[i]);
    PopHist_print_comment(chain->state->ph[0], "  0:", fp);
    PopHist_print_comment(chain->state->ph[1], "  1:", fp);
}

/**
 * Return values corresponding to the best fit.
 *
 * @param[out] bestLnObj Pointer to double, into which the function
 * will write the value of the log objective function.
 *
 * @param[out] ph Pointer to a PopHist, into which the optimal
 * population history will be written.
 */
void Chain_bestFit(Chain * chain, double *bestLnObj, PopHist * bestPh) {
    *bestLnObj = chain->bestLnObj;
    PopHist_copy(bestPh, chain->bestPh);
}

void       *runChain(void *arg) {
    Chain      *chain = (Chain *) arg;

    double      tryLnObj, mh;
    unsigned    rep, itr;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

    gsl_rng_set(rng, (unsigned) time(NULL));

    for(rep = 0; rep < chain->nreps; rep += chain->swapInterval) {

        /*
         * Run several iterations of chain (the number controlled by
         * chain->swapInterval) before attempting to swap with next
         * warmer chain.
         */
        for(itr = 0; itr < chain->swapInterval; ++itr) {
            PopHist_perturb(chain->state->ph[1], chain->state->ph[0],
                            chain->dt, chain->dNinv, rng);
            tryLnObj = lnObjFun(chain->state->ph[1], chain->u, chain->ode,
                                chain->nbins, chain->sigdsq, chain->c);

            mh = tryLnObj - chain->state->lnObj[0];

            /*
             * If flattness <1, raise mh to power "flatness",
             * which flattens the surface and facilitates movement from
             * one peak to another.
             */
            assert(chain->flatness > 0.0);
            if(chain->flatness < 1.0)
                mh *= chain->flatness;

            /* Metropolis */
            if(mh > 0.0 || log(gsl_rng_uniform_pos(rng)) <= mh) {
                PopHist    *hold = chain->state->ph[0];

                chain->state->ph[0] = chain->state->ph[1];
                chain->state->ph[1] = hold;
                chain->state->accept = 1;
                chain->state->lnObj[0] = tryLnObj;
                ++chain->naccpt;
            } else
                chain->state->accept = 0;

            if(chain->state->accept && tryLnObj > chain->bestLnObj) {
                chain->bestLnObj = tryLnObj;
                PopHist_copy(chain->bestPh, chain->state->ph[1]);
            }
        }

        /* Try swapping states with next warmer chain */
        if(chain->next) {
            mh = Chain_flatness(chain->next) - Chain_flatness(chain);

            /*
             * We want a ratio of form
             *
             *   ratio = (x^b * y^a) / (x^a * y^b)
             *
             * where x and y are the objective function values of
             * chains i and i+1, "a" is Flatness[i], and "b" is
             * Flatness[i+1]. This reduces to
             *
             *   ratio = (x/y)^(b-a)
             *
             * On log scale, this becomes
             *
             *   mh = (b-a)*(log x - log y)
             *
             */

            /* Wait on next warmer chain; LOCK WARMER CHAIN */
            PRSTAT("waiting for data from warmer chain");
            Chain_waitForData(chain->next);
            PRSTAT("got data");

            /*
             * Before comparing states on adjacent chains, we need locks
             * on both of them. Otherwise, own chain might swap state
             * with next cooler chain in the middle of the process.
             */
            Chain_lock(chain);
            mh *= chain->state->lnObj[0] - Chain_lnObj(chain->next);
            if(mh > 0.0 || log(gsl_rng_uniform_pos(rng)) < mh)
                Chain_swapState(chain);
            Chain_unlock(chain);

            Chain_setDataNeeded(chain->next);
            Chain_unlock(chain->next);  /* UNLOCK WARMER CHAIN */
            Chain_signal(chain->next);  /* wake up warmer chain */
        }

        /* Wait on next cooler chain. LOCK OWN CHAIN */
        PRSTAT("waiting until data needed");
        Chain_waitUntilDataNeeded(chain);
        PRSTAT("data needed");

        Chain_publish(chain);

        Chain_unlock(chain);    /* UNLOCK OWN CHAIN */
        Chain_signal(chain);
    }

    PRSTAT("done");

    gsl_rng_free(rng);
    pthread_exit(NULL);
    return 0;
}

/**
 * Printer header for markov chain output lines
 */
void Chain_printHdr(Chain * chain, FILE * fp) {
    int         i, n = PopHist_nParams(chain->state->ph[0]);
    char        buff[20];

    Chain_lockStdout();
    fprintf(fp, "%10s ", "lnObjFun");
    for(i = 0; i < n; ++i) {
        PopHist_paramName(chain->state->ph[0], buff, sizeof(buff), i);
        fprintf(fp, "%11s", buff);
        if(i + 1 < n)
            putc(' ', fp);
    }
    putc('\n', fp);
    Chain_unlockStdout();
}

#if 0
void Chain_printHdr(Chain * chain, FILE * fp) {
    int         i;
    char        buff1[20], buff2[20];

    Chain_lockStdout();
    fprintf(fp, "%10s ", "lnObjFun");
    for(i = 0; i < -1 + PopHist_nepoch(chain->state->ph[0]); ++i) {
        snprintf(buff1, sizeof(buff1), "1/2N[%d]", i);
        snprintf(buff2, sizeof(buff2), "duration[%d]", i);
        fprintf(fp, "%10s %11s ", buff1, buff2);
    }
    snprintf(buff1, sizeof(buff1), "1/2N[%d]", i);
    fprintf(fp, "%10s\n", buff1);
    Chain_unlockStdout();
}
#endif

/** Printer state of markov chain */
void Chain_printState(Chain * chain, FILE * fp) {
    int         i, n = PopHist_nParams(chain->state->ph[0]);

    fprintf(fp, "%10.5lg ", Chain_lnObj(chain));
    for(i = 0; i < n; ++i) {
        fprintf(fp, "%11.4lg", PopHist_paramValue(chain->state->ph[0], i));
        if(i + 1 < n)
            putc(' ', fp);
    }
    putc('\n', fp);
}
