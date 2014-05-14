/* 
 * Copyright (C) 2014 Alan Rogers
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include "annealsched.h"

/**
 * This structure represents a geometric schedule of annealing
 * temperatures.
 */
struct AnnealSched {
    int         nT;                /* number of temperature values */
    int         iT;                /* index of current tmptr */
	double      T;                 /* current temperature */
    double      initT;             /* initial temperature */

	/*
	 * temperature decays as t[i+1] = t[i]*alpha - beta
	 * t[0] is initT; t[nT-1] is 0.
	 */
    double      alpha, beta;
};

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nT      : the number of temperatures
 * initT: initial relative temperature
 * alpha   : controls rate of temperature decay. Must be in [0, 1].
 *
 * For values of alpha near 1:
 *
 * beta = 1/(n-1) + n*(a-1)/(2*(n-1)) + (a-1)^2 * n * (n-2)/(12*(n-1));
 */
AnnealSched *AnnealSched_alloc(int nT, double initT, double alpha) {
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    if(nT == 1)
        initT = 0.0;
    s->initT = initT;
    s->nT = nT;
    s->alpha = alpha;
    if(alpha <= 0.0) {
        fprintf(stderr, "%s:%d:%s alpha must be > 0\n", 
			__FILE__, __LINE__, __func__);
        exit(GSL_EINVAL);
    }else if(alpha > 1.0) {
        fprintf(stderr, "%s:%d:%s alpha must be <= 1\n", 
			__FILE__, __LINE__, __func__);
        exit(GSL_EINVAL);
    }else if(alpha == 1.0) {
        /* limit as alpha -> 1 */
        s->beta = 1.0/(nT-1.0);
    }else if(alpha < 0.99995){
        /* alpha far from 1: use full formula */
        s->beta = (alpha-1.0)/(1.0 - pow(alpha, 1.0-nT));
    }else {
        /* alpha close to 1: use 2nd-order Taylor approximation */
        double am1 = alpha - 1.0;
        s->beta = 1.0/(nT-1)
            + nT*am1/(2*(nT-1))
            + am1*am1 * (nT * (nT-2))/(12.0*(nT-1));
    }
    s->beta *= initT;
	AnnealSched_reset(s);
	return s;
}

/** Return number of temperature values */
int AnnealSched_size(const AnnealSched *sched) {
    return sched->nT;
}

void AnnealSched_reset(AnnealSched *s) {
    s->iT = 0;             /* range: 0..(nT-1) */
    s->T = s->initT;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s);
}

AnnealSched * AnnealSched_copy(const AnnealSched *old) {
    AnnealSched *new = malloc(sizeof(*new));
    if(new == NULL) {
        fprintf(stderr,"%s%d%s: bad malloc\n",
                __FILE__,__LINE__,__func__);
        exit(ENOMEM);
    }
    assert(sizeof(*old) == sizeof(*new));
    memcpy(new, old, sizeof(*new));
    return new;
}

int AnnealSched_cmp(const AnnealSched *s1, const AnnealSched *s2) {
    return memcmp( (const void *) s1, (const void *) s2,
                   sizeof(AnnealSched));
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s) {
    double tmptr = s->T;

    ++s->iT;
    if(s->iT < s->nT-1) 
        s->T = (s->T)*(s->alpha) - s->beta;
	else
        s->T = 0.0;

    return tmptr;
}

void AnnealSched_print(AnnealSched *s, FILE *fp) {
    fprintf(fp, "%s: iT=%d/%d",  __func__, s->iT, s->nT);
    fprintf(fp, " initT=%lf currT=%lf alpha=%lf beta=%lf\n",
            s->initT, s->T, s->alpha, s->beta);
}
