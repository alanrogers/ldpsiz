/**
 * @file chain.h
 * @author Alan R. Rogers
 * @brief Header for chain.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_CHAIN
#define LDPSIZ_CHAIN

#include "typedefs.h"
#include <stdio.h>

/* Turn this on for detailed status messages. */
#if 0
#define PRSTAT(MSG) \
    Chain_lockStdout();\
    printf("%s@%s:%d: %s\n",__func__,__FILE__,__LINE__, (MSG));\
    fflush(stdout); \
    Chain_unlockStdout()
#else
#define PRSTAT(MSG)
#endif

Chain      *Chain_new(unsigned which, unsigned nChains, unsigned nreps,
                      double u, double dt, double dN, double temp,
                      unsigned nbins, double *sigdsq, double *c, PopHist * ph,
                      Model * model, double odeAbsTol, double odeRelTol,
                      unsigned swapInterval);
void        Chain_free(Chain * chain);
void       *runChain(void *arg);
void        Chain_waitForData(Chain * chain);
void        Chain_waitUntilDataNeeded(Chain * chain);
void        Chain_setDataAvailable(Chain * chain);
void        Chain_setDataNeeded(Chain * chain);
void        Chain_lock(Chain * chain);
void        Chain_unlock(Chain * chain);
void        Chain_lockStdout(void);
void        Chain_unlockStdout(void);
void        Chain_signal(Chain * chain);
double      Chain_flatness(Chain * chain);
double      Chain_lnObj(Chain * chain);
void        Chain_printFull(Chain * chain, FILE * fp);
void        Chain_printHdr(Chain * chain, FILE * fp);
void        Chain_printState(Chain * chain, FILE * fp);
int         Chain_accepted(const Chain * c);
void        Chain_swapState(Chain * c);
Chain      *Chain_next(Chain * chain);
unsigned    Chain_nswapped(Chain * chain);
unsigned    Chain_naccepted(Chain * chain);
unsigned    Chain_which(Chain * chain);
void        Chain_bestFit(Chain * chain, double *bestLnObj, PopHist * bestPh);
void        Chain_sanityCheck(Chain * chain, const char *file, int line);
void        Chain_prStateAddr(const Chain * chain, int ndx, const char *file,
                              int line);
double      lnObjFun(PopHist * ph, double u, ODE * ode, int nbins,
                     double *sigdsq, double *c);

#endif
