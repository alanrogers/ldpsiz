/**
 * @file hayes.c
 * @author Alan R. Rogers
 * @brief Code for Sved's (1971) formula, as used by Hayes et
 * al. 2003. Genome Research 13:635-643.  
 */
/*
 * Internet Systems Consortium License
 * 
 * Copyright (c) 2014, Alan R. Rogers <rogers@anthro.utah.edu>
 * 
 * Permission to use, copy, modify, and/or distribute this software for
 * any purpose with or without fee is hereby granted, provided that the
 * above copyright notice and this permission notice appear in all
 * copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 * AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
 * TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <assert.h>
#include "pophist.h"
#include "hayes.h"
#include "misc.h"
#include "model.h"

static inline double svedRsq(double c, double twoN, int twoNsmp);

double      Hayes_rsq_dummy(double *y, int twoNsmp);
int         Hayes_rsqEq_dummy(double *x, double twoN, double c, double u);

double Hayes_rsq_dummy(double *y, int twoNsmp) {
    eprintf("ERR@%s:%d: should not be called", __FILE__, __LINE__);
    return *y + twoNsmp;
}

int Hayes_rsqEq_dummy(double *x, double twoN, double c, double u) {
    eprintf("ERR@%s:%d: should not be called", __FILE__, __LINE__);
    return (int) *x + twoN + c + u;
}

Model      *Model_allocHayes(int twoNsmp) {
    Model      *m = malloc(sizeof(Model));

    checkmem(m, __FILE__, __LINE__);

    m->twoNsmp = twoNsmp;
    snprintf(m->lbl, sizeof(m->lbl), "Hayes");

    m->ld = Hayes_rsq;
    m->exactLD = Hayes_rsqNoODE;
    m->ldEq = Hayes_rsqEq;
    m->stateDim = Hayes_stateDim;
    m->stateLbl = Hayes_stateLbl;
    m->stateVal = HayesData_stateVal;
    m->newState = HayesData_new;
    m->freeState = HayesData_free;

    return m;
}

/**
 * Sved's (1971; also 2009, p 184) formula for E[r^2], modified to
 * account for sampling bias as suggested by Hudson (1985, Eqn 6, p 624).
 * Hudson uses "n" throughout his paper without clearly defining
 * it. On p 631, he says that "with probability 1/n the same gamete
 * will be drawn from the sample twice" (in random sampling with
 * replacement). Thus, Hudson's "n" is the number of gametes in the
 * sample, which equals twoNsmp in the notation below.
 *
 * The HAYES_MUTATION_ADJUSTMENT is based on Tenesa et al (2007,
 * Genome Research, 17(4):520-526). On p 521, these authors say:
 *
 *    For autosomal loci, Hill (1975, Theor. Pop. Biol. 8(2):117-126) 
 *    showed that, in the presence of mutation, $E(r^2) = 
 *    (10+\rho)/(22 + 13\rho + \rho^2)$, with $\rho = 4N_e c$. Since 
 *    $(22 + 13\rho + \rho^2)$ factors into $(11+\rho)(2+\rho)$, a 
 *    further approximation is $E(r^2) = 1/(2 + \rho) = 1/(2 + 4N_e
 *    c)$.
 * 
 * Hill's formula (on his p 124), however, is for sigma_d^2 rather
 * than E[r^2], so this adjustment is relevant only to the extent that
 * sigma_d^2 approximates the expectation of r^2. Tenesa's
 * approximation also takes 10 as approximately equal to 11--a fairly
 * crude standard of approximation.
 */
static inline double svedRsq(double c, double twoN, int twoNsmp) {

#ifdef HAYES_MUTATION_ADJUSTMENT
    /* accounts for mutation */
    double      rsq = 1.0 / (2.0 + 2.0 * twoN * c);
#else
    /* no mutation */
    double      rsq = 1.0 / (1.0 + 2.0 * twoN * c);
#endif

    if(twoNsmp > 0)
        rsq += 1.0 / twoNsmp;

#if 0
    fprintf(stderr, "%s(c=%lg, twoN=%lg, twoNsmp=%d) = %lg\n",
            __func__, c, twoN, twoNsmp, rsq);
#endif

    return rsq;
}

double Hayes_rsqNoODE(double c, double u, PopHist *ph, int twoNsmp,
                      void *notused) {
    return Hayes_rsq(NULL, c, u, ph, twoNsmp);
}

double Hayes_rsq(ODE * notused, double c, double u, PopHist * ph, int twoNsmp) {
    int         j, nepochs;
    double      twoN;
    double      duration;

    myassert(ph);
    nepochs = PopHist_nepoch(ph);
    myassert(nepochs > 0);

    double      t = 1.0 / (2.0 * c);

    /* find population size at time t */
    double      thusfar = 0.0;

    for(j = 0; j < nepochs; ++j) {
        duration = PopHist_duration(ph, j);

        /* Always true on last epoch, because duration is inf */
        if(thusfar + duration > t)
            break;

        thusfar += duration;
    }
    assert(j < nepochs);
    twoN = PopHist_twoN(ph, j);

#if 0
    fprintf(stderr,"%s(c=%g, u=%g, twoNsmp=%d): t=%lg, ep=%d, twoN=%lg\n",
            __func__, c, u, twoNsmp, t, j, twoN);
#endif

    return svedRsq(c, twoN, twoNsmp);
}

double Hayes_rsqEq(double c, double u, PopHist * ph, unsigned whichEpoch,
                   int twoNsmp, void *notused) {
    double      twoN = PopHist_twoN(ph, whichEpoch);

    return svedRsq(c, twoN, twoNsmp);
}

size_t Hayes_stateDim(void) {
    return 0;
}

const char *Hayes_stateLbl(unsigned notused) {
    eprintf("ERR@%s:%d: Hayes_stateLbl should never be called",
            __FILE__, __LINE__);
    return NULL;
}

void       *HayesData_new(void) {
    return NULL;
}

double HayesData_stateVal(void *notused, unsigned i) {
    eprintf("ERR@%s:%d: HayesData_state should never be called",
            __FILE__, __LINE__);
    return 0.0;
}

void HayesData_free(void *notused) {
    myassert(notused == NULL);
    /* noop */
    return;
}
