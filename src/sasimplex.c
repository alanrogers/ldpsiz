/* multimin/sasimplex.c
 *
 * Copyright (C) 2014 Alan Rogers
 * Copyright (C) 2007, 2008, 2009 Brian Gough
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
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

/*
 * 2014-04-07: sasimplex.c and sasimplex.h implement Simplex Simulated
 *             Annealing, as described in Numerical Recipes, by Press
 *             et al. This implementation does not use code from
 *             Numerical Recipes. It is based on simplex2.c in
 *             version 1.16 of the Gnu Scientific Library, rewritten
 *             to use the adaptive simplex algorithm described by
 *             Fuchang Gao and Lixing Han. 2012. "Implementing the
 *             Nelder-Mead simplex algorithm with adaptive
 *             parameters",  (Computational Optimization and
 *             Applications  51(1):259-277, 2012).
 *
 *             Alan R. Rogers <rogers@anthro.utah.edu>
 *
 ******************************************************************
 * Documentation from simplex2.c:
 * - Originally written by Tuomo Keskitalo <tuomo.keskitalo@iki.fi>
 * - Corrections to nmsimplex_iterate and other functions
 *   by Ivo Alxneit <ivo.alxneit@psi.ch>
 * - Additional help by Brian Gough <bjg@network-theory.co.uk>
 * - Optimisations added by Brian Gough <bjg@network-theory.co.uk>
 *       + use BLAS for frequently-called functions
 *       + keep track of the center to avoid unnecessary computation
 *       + compute size as RMS value, allowing linear update on each step
 *         instead of recomputing from all N+1 vectors.
 *
 * The Simplex method of Nelder and Mead, also known as the polytope
 * search alogorithm.  Ref: Nelder, J.A., Mead, R., Computer Journal 7
 *  (1965) pp. 308-313.
 *
 * This implementation uses n+1 corner points in the simplex.
 */

#include "sasimplex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_double.h>
#include <execinfo.h>

#if 0
#define DEBUG
#else
#undef DEBUG
#endif

#define DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

/* Abort if random number seed is not yet set */
#ifndef NDEBUG
#  define ASSERT_SEED_SET(s) do{ if(!((s)->seedSet)) { \
    fprintf(stderr, "\nERR@%s:%d: call sasimplex_random_seed before %s\n", \
            __FILE__, __LINE__, __func__);                                 \
    exit(GSL_EINVAL); } }while(0)
#else
#  define ASSERT_SEED_SET(s)
#endif

#define REQUIRE(x,file,lineno,func) do {            \
  if (!(x)) { \
      dostacktrace(__FILE__,__LINE__,__func__,stderr);          \
      fprintf(stderr,"ERR@%s:%d:%s->%s:%d: Sanity check FAIL\n",\
              (file),(lineno),(func),__FILE__,__LINE__);        \
      exit(GSL_ESANITY);                                        \
   }\
} while(0)

typedef struct sasimplex_state_t sasimplex_state_t;

static inline double ran_uni(unsigned *seed);
static inline double ran_expn(unsigned *seed, double mean);
static double trial_point(const double p,
                          gsl_vector * trial,
                          const gsl_vector * h,
                          const gsl_vector * c,
                          const gsl_vector * lbound,
                          const gsl_vector * ubound,
                          const gsl_multimin_function * f);
static void update_point(sasimplex_state_t * state, size_t i,
                         const gsl_vector * x, double val);
static int  contract_by_best(sasimplex_state_t * state, double delta,
                             size_t best, gsl_vector * xc,
                             gsl_multimin_function * f);
static int  compute_center(const sasimplex_state_t * state,
                           gsl_vector * center);
static double compute_size(sasimplex_state_t * state,
                           const gsl_vector * center);
static int  sasimplex_alloc(void *vstate, size_t n);
static void sasimplex_free(void *vstate);
static int  sasimplex_set(void *vstate, gsl_multimin_function * f,
                          const gsl_vector * x,
                          double *size, const gsl_vector * step_size);
static int  sasimplex_onestep(void *vstate, gsl_multimin_function * f,
                              gsl_vector * x, double *size, double *fval);
static size_t vector_min_index(const gsl_vector * v);
static inline double sasimplex_size(sasimplex_state_t * state);
static inline double constrain_value(double x, double lo, double hi);
static void constrain_vector(gsl_vector * v, const gsl_vector * lbound,
                             const gsl_vector * ubound);
int         constrain_simplex(gsl_matrix * x1,
                              gsl_vector * fvals,
                              gsl_vector * lbound,
                              gsl_vector * ubound, gsl_multimin_function * f);
static inline int expand_around_best(sasimplex_state_t * state,
                                     gsl_multimin_function * func);
static void sasimplex_sanityCheck(const sasimplex_state_t * state,
                                  const char *file, int lineno,
                                  const char *func);
static void dostacktrace(const char *file, int line, const char *func,
                         FILE * ofp);
static void vector_print(gsl_vector *x, FILE *fp);
static void sasimplex_print(sasimplex_state_t * state);
static const char *statusLbl(int status);


/** The state of the minimizer */
struct sasimplex_state_t {
    gsl_matrix *x1;             /* simplex corner points */
    gsl_vector *f1;             /* function values at corner points */
    gsl_vector *ws1;            /* workspace 1 for algorithm */
    gsl_vector *ws2;            /* workspace 2 for algorithm */
    gsl_vector *center;         /* center of all points */
    gsl_vector *delta;          /* current step */
    gsl_vector *xmc;            /* x - center (workspace) */
    gsl_vector *lbound;         /* lower bound */
    gsl_vector *ubound;         /* upper bound */
    gsl_vector *step_size;      /* size of initial simplex */
    double      S2;
    double      temperature;    /* increase to flatten surface */
    double      bestEver;       /* best func val ever seen */
    unsigned    seedSet;        /* 0 initially; 1 after seed is set */
    unsigned    seed;           /* for random number generator */
};

static const gsl_multimin_fminimizer_type sasimplex_type = {
    "sasimplex",                /* name */
    sizeof(sasimplex_state_t),
    &sasimplex_alloc,
    &sasimplex_set,
    &sasimplex_onestep,
    &sasimplex_free
};

const       gsl_multimin_fminimizer_type
    * gsl_multimin_fminimizer_sasimplex = &sasimplex_type;

static const char *statusLbl(int status) {
            switch(status) {
            case GSL_SUCCESS:
                return("success");
            case GSL_ETOLX:
                return("etolx");
            case GSL_ETOLF:
                return("etolf");
            case GSL_CONTINUE:
                return("continue");
            default:
                return("????");
            }
}


/** Print an array of doubles to file fp. len is length of array */
static void vector_print(gsl_vector *x, FILE *fp) {
    size_t i;

    fputs("[", fp);
    for(i=0; i < x->size; ++i) {
        fprintf(fp, "%lf", gsl_vector_get(x, i));
        if(i+1 < x->size)
            fputs(", ", fp);
    }
    fputs("]", fp);
}

static void dostacktrace(const char *file, int line, const char *func,
                         FILE * ofp) {
    enum { CALLSTACK_SIZE = 128 };
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace returned %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:%s\n", file, line, func);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
}

static void sasimplex_sanityCheck(const sasimplex_state_t * state,
                                  const char *file, int lineno,
                                  const char *func) {
#ifndef NDEBUG
    REQUIRE(state->x1 != NULL, file, lineno, func);
    REQUIRE(state->f1 != NULL, file, lineno, func);
    REQUIRE(state->ws1 != NULL, file, lineno, func);
    REQUIRE(state->ws2 != NULL, file, lineno, func);
    REQUIRE(state->center != NULL, file, lineno, func);
    REQUIRE(state->delta != NULL, file, lineno, func);
    REQUIRE(state->center != NULL, file, lineno, func);

    /* either both NULL or both non-NULL */
    REQUIRE((state->lbound == NULL) == (state->ubound == NULL),
            file, lineno, func);

    size_t      i, j;
    size_t      n = state->center->size;
    REQUIRE(n > 0, file, lineno, func);

    if(state->lbound != NULL) {
        double      lb, ub, x;
        /* Is simplex within bounds? */
        for(i = 0; i <= n; ++i) {
            gsl_vector_const_view row = gsl_matrix_const_row(state->x1, i);
            for(j = 0; j < n; ++j) {
                lb = gsl_vector_get(state->lbound, j);
                ub = gsl_vector_get(state->ubound, j);
                x = gsl_vector_get(&row.vector, j);
                REQUIRE(lb <= x, file, lineno, func);
                REQUIRE(x <= ub, file, lineno, func);
            }
        }

        /* Is center within bounds? */
        for(j = 0; j < n; ++j) {
            lb = gsl_vector_get(state->lbound, j);
            ub = gsl_vector_get(state->ubound, j);
            lb -= 2.0 * DBL_EPSILON;
            ub += 2.0 * DBL_EPSILON;
            x = gsl_vector_get(state->center, j);
            if(x < lb) {
                printf("%s:%d:%s: x=%lf < %lf; diff=%le]\n",
                       __FILE__, __LINE__, __func__, x, lb, lb - x);
            } else if(x > ub) {
                printf("%s:%d:%s: x=%lf > %lf; diff=%le]\n",
                       __FILE__, __LINE__, __func__, x, ub, x - ub);
            }
            REQUIRE(lb <= x, file, lineno, func);
            REQUIRE(x <= ub, file, lineno, func);
        }
    }

    /* Is center really the center? */
    (void) compute_center(state, state->ws1);
    double      absdiff = 0.0, absmax = 0.0;
    const double reltol = 0.01;
    for(j = 0; j < n; ++j) {
        double      x = gsl_vector_get(state->ws1, j);
        double      y = gsl_vector_get(state->center, j);
        absdiff += fabs(x - y);
        absmax += fmax(fabs(x), fabs(y));
    }
    if(absdiff > absmax * reltol) {
        printf("absdiff=%le absmax=%le relerr=%le > %le\n",
               absdiff, absmax, absdiff / absmax, reltol);
        fflush(stdout);
    }
    REQUIRE(absdiff <= absmax * reltol, file, lineno, func);
#endif
}

/** Print state variable of fminimizer */
void fminimizer_print(gsl_multimin_fminimizer * minimizer) {
    sasimplex_state_t *state = minimizer->state;
    sasimplex_print(state);
}

/** Print object of type sasimplex_state_t */
static void sasimplex_print(sasimplex_state_t *state) {
    size_t      n = state->center->size;    /* dimension of state vector */
    size_t      i, j;

    printf("Simplex\n");
    for(i = 0; i <= n; ++i) {
        for(j = 0; j < n; ++j)
            printf(" %lf", gsl_matrix_get(state->x1, i, j));
        printf(": %lf\n", gsl_vector_get(state->f1, i));
    }
    printf(" tmptr=%lf bestEver=%lf\n", state->temperature, state->bestEver);
#if 0
    if(state->lbound != NULL) {
        printf(" upper bound:");
        for(i = 0; i < n; ++i)
            printf(" %lf", gsl_vector_get(state->ubound, i));
        putchar('\n');
    }
    if(state->ubound != NULL) {
        printf(" lower bound:");
        for(i = 0; i < n; ++i)
            printf(" %lf", gsl_vector_get(state->lbound, i));
        putchar('\n');
    }
#endif
}

/*
 * This code is based on the source for gsl_vector_min_index, which
 * returns the index of the smallest entry in vector v.
 *
 * It behaves exactly like gsl_vector_min_index, provided that the
 * latter is compiled with macro FP undefined: it returns the index of
 * the minimum non-NaN value, ignoring NaN entries. If all entries are
 * NaN, or if the vector has zero length, it returns 0.
 *
 * The present code is an effort to avoid the behavior that
 * gsl_min_index exhibits if compiled with macro FP defined. In that
 * case, if any NaN values are present in the vector, gsl_min_index
 * returns the index of the 1st NaN it encounters. This causes the
 * simplex algorithm to gravitate toward regions of the state space
 * where the function value is undefined.
 *
 * The present code uses gsl_vector_get to retrieve values and changes
 * the lower bound of the loop from 0 to 1.
 */
static size_t vector_min_index(const gsl_vector * v) {
    const size_t N = v->size;

    double      min = gsl_vector_get(v, 0);
    size_t      imin = 0;
    size_t      i;

    for(i = 1; i < N; i++) {
        double      x = gsl_vector_get(v, i);
        if(x < min) {
            min = x;
            imin = i;
        }
    }

    return imin;
}

/** Set temperature. */
void
sasimplex_set_temp(gsl_multimin_fminimizer * minimizer, double temperature) {
    sasimplex_state_t *state = minimizer->state;

    state->temperature = temperature;
}

/** Test convergence based on relative spread of function values */
int
sasimplex_converged(gsl_multimin_fminimizer * minimizer, double *size,
                    double tol_fval, double tol_size) {
    sasimplex_state_t *state = minimizer->state;

    double      best, worst, dy, tol1;
    int         status = GSL_ESANITY;
    int         fvals_eq;    

#if 0
    gsl_vector_minmax(state->f1, &best, &worst);
#else
    worst = gsl_vector_max(state->f1);
    best = state->bestEver;
#endif
    assert(worst >= best);
    /* convergence criteria from zeroin */
    tol1 = 4.0 * tol_fval * (fabs(worst) + fabs(best)) + tol_fval;
    dy = fabs(worst - best);

    fvals_eq = (dy < tol1 ? 1 : 0);
    *size = sasimplex_size(state);
    if(*size < tol_size) {
        if(fvals_eq) {
            status = GSL_SUCCESS;
        } else {
            /* stuck */
            status = GSL_ETOLF;
#ifdef DEBUG
            sasimplex_print(state);
#endif
            DPRINTF(("%s:%d:%s: stuck. worst=%lf best=%lf diff=%lf >= %lf = tol1\n",
                     __FILE__,__LINE__,__func__,
                     worst, best, dy, tol1)); 
        }
    } else {
        if(fvals_eq) {
            /* flat objective function */
            status = GSL_ETOLX;
        } else {
            /* keep going */
            status = GSL_CONTINUE;
        }
    }

    assert(status != GSL_ESANITY);

    return status;
}

/**
 * Fold x back and forth across the boundaries "lo" and "hi" to obtain a value
 * y such that lo <= y <= hi.
 */
static inline double constrain_value(double x, double lo, double hi) {
    assert(hi > lo);
    x -= lo;
    hi -= lo;
    double      z = fabs(fmod(x, 2.0 * hi));
    if(z > hi)
        z = 2.0 * hi - z;
    z += lo;
    assert(z >= lo && z <= hi + lo);
    return z;
}

/** Contrain all entries of v to lie between lower and upper bounds */
static void constrain_vector(gsl_vector * v,
                             const gsl_vector * lbound,
                             const gsl_vector * ubound) {
    size_t      i;
    assert((lbound && ubound) || ((!lbound) && (!ubound)));
    if(lbound != NULL) {
        for(i = 0; i < v->size; ++i) {
            double      val = gsl_vector_get(v, i);
            double      lb = gsl_vector_get(lbound, i);
            double      ub = gsl_vector_get(ubound, i);
            val = constrain_value(val, lb, ub);
            gsl_vector_set(v, i, val);
        }
    }
}

/*
 * Calculate a new trial vertex and associated function value.
 *
 * On input:
 *
 * p : the coefficient that determines whether we're doing a reflection,
 *     an extension, an outside contraction, or an inside contraction.
 * h : the vertex of simplex with highest function value
 * c : the centroid of the simplex, including all n+1 points
 *
 * On output:
 * trial :  the new trial vectex
 *
 * Function returns the function value at the new vertex.
 */
static double trial_point(const double p, gsl_vector * trial,   /* to hold new vertex */
                          const gsl_vector * h, /* highest point in simplex */
                          const gsl_vector * c, /* centroid of entire simplex */
                          const gsl_vector * lbound,    /* lower bound */
                          const gsl_vector * ubound,    /* lower bound */
                          const gsl_multimin_function * f) {

    /* n is the dimension of the state space */
    const size_t n = h->size;

    /*
     * We want to calculate a trial vertex as
     *
     *    trial = (1-p)*m + p*h = m - p*(m - h)
     *
     * where h is the vector defining the simplex point with
     * highest function value and m is the centroid of the
     * remaining points. However, the present code calculates c,
     * the center of the entire simplex, not excluding the highest
     * point. The formula above for "trial" is equivalent to
     *
     *    trial = (1-a)*c + a*h = c - a*(c - h)
     *
     * where
     *
     *      a = ((n+1)*p - 1)/n
     */
    double      a = ((n + 1) * p - 1.0) / n;
    size_t      i;
    int         ineq_constraint = 0;

    assert((lbound && ubound) || ((!lbound) && (!ubound)));

    if(lbound != NULL)
        ineq_constraint = 1;

    /* trial = (1-a)*c + a*h */
    for(i = 0; i < n; ++i) {
        double      cval = gsl_vector_get(c, i);
        double      hval = gsl_vector_get(h, i);
        double      try = cval - a * (cval - hval);
        if(ineq_constraint) {
            double      lb = gsl_vector_get(lbound, i);
            double      ub = gsl_vector_get(ubound, i);
            try = constrain_value(try, lb, ub);
            assert(try >= lb);
            assert(try <= ub);
        }
        gsl_vector_set(trial, i, try);
    }
#if 0
    gsl_vector_memcpy(trial, c);
    gsl_blas_dscal(1.0 - a, trial);
    gsl_blas_daxpy(a, h, trial);
#endif

    /* new function value */
    double      newval = GSL_MULTIMIN_FN_EVAL(f, trial);

    return newval;
}

static void
update_point(sasimplex_state_t * state, size_t i,
             const gsl_vector * x, double val) {
    gsl_vector_const_view x_orig = gsl_matrix_const_row(state->x1, i);
    const size_t P = state->x1->size1;

    /* Compute state->delta = x - x_orig */
    gsl_vector_memcpy(state->delta, x);
    gsl_blas_daxpy(-1.0, &x_orig.vector, state->delta);

    /* Compute state->xmc = x_orig - c */
    gsl_vector_memcpy(state->xmc, &x_orig.vector);
    gsl_blas_daxpy(-1.0, state->center, state->xmc);

    /* Update size: S2' = S2 + (2/P) * (x_orig - c).delta + (P-1)*(delta/P)^2 */
    {
        double      d = gsl_blas_dnrm2(state->delta);
        double      xmcd;

        /* increments state->S2 */
        gsl_blas_ddot(state->xmc, state->delta, &xmcd);
        state->S2 += (2.0 / P) * xmcd + ((P - 1.0) / P) * (d * d / P);
    }

    /* Update state->center:  c' = c + (x - x_orig) / P */
    {
        double      alpha = 1.0 / P;

        /* result goes in state->center */
        gsl_blas_daxpy(-alpha, &x_orig.vector, state->center);
        gsl_blas_daxpy(alpha, x, state->center);
    }

    gsl_matrix_set_row(state->x1, i, x);
    gsl_vector_set(state->f1, i, val);
    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);
}

static inline int
expand_around_best(sasimplex_state_t * state, gsl_multimin_function * func) {
    DPRINTF(("%s\n", __func__));
    size_t      best = vector_min_index(state->f1);
    double      inflateBy = 200.0;
    int         status;

    status = contract_by_best(state, inflateBy, best, state->ws1, func);
    if(status != GSL_SUCCESS)
        GSL_ERROR("contract_by_best failed", status);
    return status;
}

/*
 * Function contracts the simplex toward current minimum. That is, all
 * corners besides the best corner are moved. (This function is rarely
 * called in practice, since it is the last choice, hence not
 * optimised - BJG)
 */
static int
contract_by_best(sasimplex_state_t * state, double delta,
                 size_t best, gsl_vector * xc, gsl_multimin_function * f) {
    DPRINTF(("%s\n", __func__));
    /* the xc vector is simply work space here */
    gsl_matrix *x1 = state->x1;
    gsl_vector *f1 = state->f1;
    size_t      i, j;
    int         status = GSL_SUCCESS;

    for(i = 0; i < x1->size1; i++) {
        if(i == best)
            continue;
        double      newval;
        for(j = 0; j < x1->size2; j++) {
            double      xbest = gsl_matrix_get(x1, best, j);
            double      xi = gsl_matrix_get(x1, i, j);
            newval = xbest + delta * (xi - xbest);
            gsl_matrix_set(x1, i, j, newval);
        }

        /* evaluate function in the new point */
        gsl_matrix_get_row(xc, x1, i);
        newval = GSL_MULTIMIN_FN_EVAL(f, xc);
        gsl_vector_set(f1, i, newval);

        if(gsl_finite(newval)) {
            if(newval < state->bestEver) {
                state->bestEver = newval;
            }
        }else
            GSL_ERROR("bad function value", GSL_EBADFUNC);
    }

    if(delta > 1.0) {
        status = constrain_simplex(state->x1, state->f1, state->lbound,
                                   state->ubound, f);
        if(status != GSL_SUCCESS)
            GSL_ERROR("constrain_simplex failed", status);
        
    }

    /* We need to update the centre and size as well */
    compute_center(state, state->center);
    compute_size(state, state->center);

    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);

    return status;
}

/*
 * Calculate the center of the simplex and stores in center.
 *
 * This routine calculates the center of all points, rather than
 * of only the P-1 best points. The difference is only cosmetic,
 * because one can express the algorithm either way. It is important
 * to keep straight, however, when trying to reconcile variants of the
 * Nelder-Mead algorithm.
 */
static int
compute_center(const sasimplex_state_t * state, gsl_vector * center) {
    gsl_matrix *x1 = state->x1;
    const size_t P = x1->size1;
    size_t      i;

    gsl_vector_set_zero(center);

    for(i = 0; i < P; i++) {
        gsl_vector_const_view row = gsl_matrix_const_row(x1, i);
        gsl_blas_daxpy(1.0, &row.vector, center);
    }
    gsl_blas_dscal(1.0 / P, center);
    return GSL_SUCCESS;
}

static inline double sasimplex_size(sasimplex_state_t * state) {
    double      S2 = state->S2;

    if(S2 > 0.0)
        return sqrt(S2);
    /* recompute if accumulated error has made size invalid */
    return compute_size(state, state->center);
}

static double
compute_size(sasimplex_state_t * state, const gsl_vector * center) {
    /* calculates simplex size as rms sum of length of vectors
     * from simplex center to corner points:
     *
     * sqrt( sum ( || y - y_middlepoint ||^2 ) / n )
     */
    gsl_vector *s = state->ws1;
    gsl_matrix *x1 = state->x1;
    const size_t P = x1->size1;
    size_t      i;
    double      ss = 0.0;

    for(i = 0; i < P; i++) {
        double      t;
        gsl_matrix_get_row(s, x1, i);
        gsl_blas_daxpy(-1.0, center, s);
        t = gsl_blas_dnrm2(s);
        ss += t * t;
    }

    /* Store squared size in the state */
    state->S2 = (ss / P);

    return sqrt(ss / P);
}

/**
 * Allocate arrays within an object of type sasimplex_state_t.
 * The object itself must be allocated previously.
 */
static int sasimplex_alloc(void *vstate, size_t n) {
    sasimplex_state_t *state = (sasimplex_state_t *) vstate;

    if(n == 0) {
        GSL_ERROR("invalid number of parameters specified", GSL_EINVAL);
    }

    state->x1 = gsl_matrix_alloc(n + 1, n);
    if(state->x1 == NULL) {
        GSL_ERROR("failed to allocate space for x1", GSL_ENOMEM);
    }

    state->f1 = gsl_vector_alloc(n + 1);
    if(state->f1 == NULL) {
        gsl_matrix_free(state->x1);
        GSL_ERROR("failed to allocate space for y", GSL_ENOMEM);
    }

    state->ws1 = gsl_vector_alloc(n);
    if(state->ws1 == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        GSL_ERROR("failed to allocate space for ws1", GSL_ENOMEM);
    }

    state->ws2 = gsl_vector_alloc(n);
    if(state->ws2 == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        gsl_vector_free(state->ws1);
        GSL_ERROR("failed to allocate space for ws2", GSL_ENOMEM);
    }

    state->center = gsl_vector_alloc(n);
    if(state->center == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        gsl_vector_free(state->ws1);
        gsl_vector_free(state->ws2);
        GSL_ERROR("failed to allocate space for center", GSL_ENOMEM);
    }

    state->delta = gsl_vector_alloc(n);
    if(state->delta == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        gsl_vector_free(state->ws1);
        gsl_vector_free(state->ws2);
        gsl_vector_free(state->center);
        GSL_ERROR("failed to allocate space for delta", GSL_ENOMEM);
    }

    state->xmc = gsl_vector_alloc(n);
    if(state->xmc == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        gsl_vector_free(state->ws1);
        gsl_vector_free(state->ws2);
        gsl_vector_free(state->center);
        gsl_vector_free(state->delta);
        GSL_ERROR("failed to allocate space for xmc", GSL_ENOMEM);
    }

    state->step_size = gsl_vector_alloc(n);
    if(state->step_size == NULL) {
        gsl_matrix_free(state->x1);
        gsl_vector_free(state->f1);
        gsl_vector_free(state->ws1);
        gsl_vector_free(state->ws2);
        gsl_vector_free(state->center);
        gsl_vector_free(state->delta);
        gsl_vector_free(state->xmc);
        GSL_ERROR("failed to allocate space for step_size", GSL_ENOMEM);
    }

    state->lbound = state->ubound = NULL;

    state->temperature = 0.0;
    state->seedSet = state->seed = 0;
    state->bestEver = DBL_MAX;

    return GSL_SUCCESS;
}

static void sasimplex_free(void *vstate) {
    sasimplex_state_t *state = (sasimplex_state_t *) vstate;

    gsl_matrix_free(state->x1);
    gsl_vector_free(state->f1);
    gsl_vector_free(state->ws1);
    gsl_vector_free(state->ws2);
    gsl_vector_free(state->center);
    gsl_vector_free(state->delta);
    gsl_vector_free(state->xmc);
    gsl_vector_free(state->step_size);
    if(state->lbound != NULL)
        gsl_vector_free(state->lbound);
    if(state->ubound != NULL)
        gsl_vector_free(state->ubound);
}

/** Set lower and upper bounds on state vector */
int sasimplex_set_bounds(gsl_multimin_fminimizer * minimizer,
                         const gsl_vector * lbound,
                         const gsl_vector * ubound,
                         const gsl_vector * step_size) {
    sasimplex_state_t *state = minimizer->state;
    int         status = GSL_SUCCESS;

    if(lbound == NULL || ubound == NULL)
        GSL_ERROR("NULL vector passed to sasimplex_set_bounds", GSL_EDOM);

    state->lbound = gsl_vector_alloc(lbound->size);
    state->ubound = gsl_vector_alloc(lbound->size);
    if(state->lbound == NULL || state->ubound == NULL)
        GSL_ERROR("bad gsl_vector_alloc", GSL_ENOMEM);

    gsl_vector_memcpy(state->lbound, lbound);
    gsl_vector_memcpy(state->ubound, ubound);

    status = constrain_simplex(state->x1, state->f1, state->lbound,
                               state->ubound, minimizer->f);
    if(status != GSL_SUCCESS)
        GSL_ERROR("constrain_simplex failed", status);

    compute_center(state, state->center);
    (void) compute_size(state, state->center);
    state->bestEver = gsl_vector_min(state->f1);

    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);
    return status;
}

/**
 * Make sure simplex points (in matrix x1) are between lbound and
 * ubound and adjust function values (in vector fvals).
 */
int constrain_simplex(gsl_matrix * x1,
                      gsl_vector * fvals,
                      gsl_vector * lbound,
                      gsl_vector * ubound, gsl_multimin_function * f) {

    if(lbound == NULL && ubound == NULL)
        return 0;

    assert(lbound && ubound);

    size_t      i, npts = fvals->size;

    for(i = 0; i < npts; ++i) {
        gsl_vector_view x = gsl_matrix_row(x1, i);
        constrain_vector(&x.vector, lbound, ubound);
        double      val = GSL_MULTIMIN_FN_EVAL(f, &x.vector);
        if(!gsl_finite(val))
            GSL_ERROR("non-finite function value encountered", GSL_EBADFUNC);
        gsl_vector_set(fvals, i, val);
    }
    return GSL_SUCCESS;
}

/**
 * Random variates drawn from standard uniform distribution
 */
static inline double ran_uni(unsigned *seed) {
    return rand_r(seed) / (RAND_MAX + 1.0);
}

/**
 * Random variates drawn from exponential distribution with given
 * mean.
 */
static inline double ran_expn(unsigned *seed, double mean) {
    double      u;
    do {
        u = rand_r(seed);
    } while(u == 0.0);
    u /= RAND_MAX;              /* u is uniform on (0,1] */
    return -mean * log(u);
}

/** Provide seed for random number generator. */
void sasimplex_random_seed(gsl_multimin_fminimizer * minimizer, unsigned seed) {
    sasimplex_state_t *state = minimizer->state;

    state->seed = seed;
    state->seedSet = 1;
}

/*
 * Measure vertical scale of simplex as the difference between the
 * current minimum function value and the smallest value ever seen.
 */
double sasimplex_vertical_scale(gsl_multimin_fminimizer * minimizer) {
    sasimplex_state_t *state = minimizer->state;
    double      currBest = gsl_vector_min(state->f1);

    assert(state->bestEver < DBL_MAX);
    assert(currBest >= state->bestEver);

#if 1
    /*
     * If currBest equals bestEver, then return the difference between
     * the current max and min function values.
     */
    if(currBest == state->bestEver)
        currBest = gsl_vector_max(state->f1);

    assert(currBest > state->bestEver);
#endif

    return currBest - state->bestEver;
}

static int
sasimplex_set(void *vstate, gsl_multimin_function * func,
              const gsl_vector * x,
              double *size, const gsl_vector * step_size) {
    double      val;
    int         status;
    sasimplex_state_t *state = (sasimplex_state_t *) vstate;
    gsl_vector *xtemp = state->ws1;

    if(xtemp->size != x->size)
        GSL_ERROR("incompatible size of x", GSL_EINVAL);

    if(xtemp->size != step_size->size)
        GSL_ERROR("incompatible size of step_size", GSL_EINVAL);

    status = gsl_vector_memcpy(state->step_size, step_size);
    if(status != 0)
        GSL_ERROR("gsl_vector_memcpy failed", GSL_EFAILED);

    /* first point is x0 */
    val = GSL_MULTIMIN_FN_EVAL(func, x);
    if(!gsl_finite(val))
        GSL_ERROR("non-finite function value", GSL_EBADFUNC);

    gsl_matrix_set_row(state->x1, 0, x);
    gsl_vector_set(state->f1, 0, val);

    size_t      i;
    gsl_vector_view x0 = gsl_matrix_row(state->x1, 0);

    /* following points are initialized to x0 + step_size */
    for(i = 0; i < x->size; i++) {
        status = gsl_vector_memcpy(xtemp, &x0.vector);
        if(status != 0)
            GSL_ERROR("gsl_vector_memcpy failed", GSL_EFAILED);

        {
            double      xi = gsl_vector_get(&x0.vector, i);
            xi += gsl_vector_get(step_size, i);
            gsl_vector_set(xtemp, i, xi);
        }
        val = GSL_MULTIMIN_FN_EVAL(func, xtemp);
        if(!gsl_finite(val))
            GSL_ERROR("non-finite function value encountered", GSL_EBADFUNC);
        gsl_matrix_set_row(state->x1, i + 1, xtemp);
        gsl_vector_set(state->f1, i + 1, val);
    }

    status = constrain_simplex(state->x1, state->f1, state->lbound,
                               state->ubound, func);
    if(status != GSL_SUCCESS)
        GSL_ERROR("constrain_simplex failed", GSL_EFAILED);
    compute_center(state, state->center);
    *size = compute_size(state, state->center);
    state->bestEver = gsl_vector_min(state->f1);

    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);
    return status;
}

/** One simplex iteration with simulated annealing. */
static int
sasimplex_onestep(void *vstate, gsl_multimin_function * func,
                  gsl_vector * x, double *size, double *fval) {

    sasimplex_state_t *state = (sasimplex_state_t *) vstate;

    /* xc and xc2 vectors store tried corner point coordinates */
    gsl_vector *xc = state->ws1;
    gsl_vector *xc2 = state->ws2;
    gsl_vector *fvec = state->f1;
    gsl_matrix *x1 = state->x1;
    const size_t n = fvec->size;
    size_t      i;
    size_t      hi, lo;
    double      dhi, ds_hi, dlo, hold;
    int         status = GSL_SUCCESS;
    double      v, v2;          /* unperturbed trial values */
    double      pv, pv2;        /* perturbed trial values */
    double      temp = state->temperature;

    if(xc->size != x->size)
        GSL_ERROR("incompatible size of x", GSL_EINVAL);

    /*
     * Constants from Eqn 4.1 of "Implementing the Nelder-Mead simplex
     * algorithm with adaptive parameters", by Fuchang Gao and Lixing
     * Han (Computational Optimization and Applications
     * 51(1):259-277, 2012).
     */
    double s = fmin(2.0, n);
    const double alpha = 1.0;
    const double beta = 1.0 + 2.0 / s;
	const double gmma = 0.75 - 1.0 / (2.0 * s);
	const double delta = 1.0 - 1.0 / s;

    /*
     * Find highest, second highest and lowest point. We need the
     * indices (lo and hi) of the  low and high points, but we don't
     * need the index of the second highest. We need the function
     * values of all three.
     *
     * dlo, ds_hi, and dhi are function values at these three points,
     * perturbed upward by random amounts. They are thus somewhat
     * worse than the true function values.
     *
     * Because of the perturbations, we can't use
     * gsl_vector_minmax_index.
     */
    lo = 0;
    hi = 1;
    ASSERT_SEED_SET(state);
    dlo = gsl_vector_get(fvec, lo) + ran_expn(&state->seed, temp);
    dhi = gsl_vector_get(fvec, hi) + ran_expn(&state->seed, temp);

    if(dhi < dlo) {             /* swap lo and hi */
        lo = 1;
        hi = 0;
        hold = dlo;
        dlo = dhi;
        dhi = hold;
    }
    ds_hi = dlo;

    for(i = 2; i < n; i++) {
        v = gsl_vector_get(fvec, i) + ran_expn(&state->seed, temp);
        if(v < dlo) {
            dlo = v;
            lo = i;
        } else if(v > dhi) {
            ds_hi = dhi;
            dhi = v;
            hi = i;
        } else if(v > ds_hi) {
            ds_hi = v;
        }
    }

    /* simplex point with worst (largest) func value */
    gsl_vector_const_view hvec = gsl_matrix_const_row(state->x1, hi);

    /*
     * v is the true function value at the trial point and pv is the
     * perturbed version of that value. In contrast to the upward
     * perturbations in dlo, ds_hi, and dhi, the perturbation here is
     * downward, making the trial value a little better from the
     * perspective of the minimizer. This encourages the algorithm to
     * accept trial values--makes it eager to explore.
     */

    /* reflection */
    v = trial_point(-alpha, xc, &hvec.vector, state->center, state->lbound,
                    state->ubound, func);
    pv = v - ran_expn(&state->seed, temp);

    if(pv < dlo) {              /* try expansion */
        v2 = trial_point(-alpha * beta, xc2, &hvec.vector, state->center,
                         state->lbound, state->ubound, func);
        pv2 = v2 - ran_expn(&state->seed, temp);
        if(pv2 < pv) {          /* accept expansion */
            update_point(state, hi, xc2, v2);
            DPRINTF(("%s:%d:%s: expansion\n", __FILE__,__LINE__,__func__));
        } else {                /* accept reflection */
            update_point(state, hi, xc, v);
            DPRINTF(("%s:%d:%s: reflection\n", __FILE__,__LINE__,__func__));
        }
    } else if(pv < ds_hi) {     /* accept reflection */
        assert(dlo <= pv);
        update_point(state, hi, xc, v);
        DPRINTF(("%s:%d:%s: reflection\n", __FILE__,__LINE__,__func__));
    } else if(pv < dhi) {       /* try outside contraction */
        assert(ds_hi <= pv);
        v2 = trial_point(-alpha * gmma, xc2, &hvec.vector, state->center,
                         state->lbound, state->ubound, func);
        pv2 = v2 - ran_expn(&state->seed, temp);
        if(pv2 <= pv) {         /* accept outside contraction */
            update_point(state, hi, xc2, v2);
            DPRINTF(("%s:%d:%s: outside contraction\n", __FILE__,__LINE__,__func__));
        } else {                /* shrink */
            status = contract_by_best(state, delta, lo, xc, func);
            if(status != GSL_SUCCESS)
                GSL_ERROR("contract_by_best failed", status);
            DPRINTF(("%s:%d:%s: shrink\n", __FILE__,__LINE__,__func__));
        }
    } else {                    /* try inside contraction */
        assert(dhi <= pv || !gsl_finite(pv));
        /*
         * According to Gao and Han (3rd page), the inside contraction is
         *
         *   center*(1-alpha*gmma)  + alpha*gmma*h
         *
         * According to Lagarias et al (eqn 2.7), it is
         *
         *   center*(1-gmma)  + gmma*h
         *
         * This shouldn't matter, because Gao and Han (eqn 4.1) take
         * alpha=1.
         */
        v2 = trial_point(alpha * gmma, xc2, &hvec.vector, state->center,
                         state->lbound, state->ubound, func);
        pv2 = v2 - ran_expn(&state->seed, temp);
        if(pv2 < dhi) {         /* accept inside contraction */
            update_point(state, hi, xc2, v2);
            DPRINTF(("%s:%d:%s: inside contraction\n", __FILE__,__LINE__,__func__));
        } else {                /* shrink */
            status = contract_by_best(state, delta, lo, xc, func);
            if(status != GSL_SUCCESS)
                GSL_ERROR("contract_by_best failed", status);
            DPRINTF(("%s:%d:%s: contract by best\n", __FILE__,__LINE__,__func__));
        }
    }

    if(status != GSL_SUCCESS)
        GSL_ERROR("contraction failed", GSL_EFAILED);

    /* Return lowest point of simplex as x. */
    lo = vector_min_index(fvec);
    gsl_matrix_get_row(x, x1, lo);
    *fval = gsl_vector_get(fvec, lo);

    if(*fval < state->bestEver) {
        state->bestEver = *fval;
    }

    *size = sasimplex_size(state);

    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);
    return status;
}

int
sasimplex_randomize_state(gsl_multimin_fminimizer * minimizer,
                          int rotate, gsl_vector * lo,
                          gsl_vector * hi, const gsl_vector * step_size) {
    sasimplex_state_t *state = minimizer->state;
    gsl_multimin_function *func = minimizer->f;
    double      val;
    size_t      i, j;
    gsl_vector *xtemp = state->ws1;
    size_t      stateDim = xtemp->size;
    int         status = GSL_SUCCESS;

    ASSERT_SEED_SET(state);
    /*
     * Copy of point 0 of the simplex into xtemp.
     *
     * It's not clear this is necessary. Current code copies row 0 into
     * xtemp, manipulates xtemp, then copies back into row 0. Perhaps
     * I could get away with working directly on pt0.
     */
    gsl_vector_const_view row0 = gsl_matrix_const_row(state->x1, 0);
    gsl_vector_memcpy(xtemp, &row0.vector);

    /*
     * If lo and hi exist, then initialize around random point.
     * Otherwise, initial point will be first row of existing
     * matrix state->x1.
     */
    if(lo != NULL && hi != NULL) {
        if(stateDim != lo->size) {
            GSL_ERROR("incompatible size of lo", GSL_EINVAL);
        }
        if(stateDim != hi->size) {
            GSL_ERROR("incompatible size of hi", GSL_EINVAL);
        }
        for(i = 0; i < stateDim; ++i) {
            double      y = gsl_vector_get(lo, i);
            double      z = gsl_vector_get(hi, i);
            assert(y <= z);
            val = y + (z - y) * ran_uni(&state->seed);
            gsl_vector_set(xtemp, i, val);
        }
        gsl_matrix_set_row(state->x1, 0, xtemp);
        gsl_vector_memcpy(minimizer->x, xtemp);
        val = GSL_MULTIMIN_FN_EVAL(func, xtemp);
        if(!gsl_finite(val)) {
            GSL_ERROR("non-finite function value encountered", GSL_EBADFUNC);
        }
        gsl_vector_set(state->f1, 0, val);
    }

    gsl_matrix_view m =
        gsl_matrix_submatrix(state->x1, 1, 0, stateDim, stateDim);

    gsl_matrix_set_identity(&m.matrix);

    /* start with random reflections */
    for(i = 0; i < stateDim; i++) {
        if(0.5 < ran_uni(&state->seed))
            gsl_matrix_set(&m.matrix, i, i, -1.0);
    }

    if(rotate) {
        /* apply random rotations */
        for(i = 0; i < stateDim; i++) {
            for(j = i + 1; j < stateDim; j++) {
                /* rotate columns i and j by a random angle */
                double      angle = 2.0 * M_PI * ran_uni(&state->seed);
                double      c = cos(angle), s = sin(angle);
                gsl_vector_view c_i = gsl_matrix_column(&m.matrix, i);
                gsl_vector_view c_j = gsl_matrix_column(&m.matrix, j);

                gsl_blas_drot(&c_i.vector, &c_j.vector, c, s);
            }
        }
    }

    /* scale the orthonormal basis by the user-supplied step_size in
     * each dimension, and use as an offset from the central point x */
    for(i = 0; i < stateDim; i++) {
        double      x_i = gsl_vector_get(&row0.vector, i);
        double      s_i = gsl_vector_get(step_size, i);
        gsl_vector_view c_i = gsl_matrix_column(&m.matrix, i);

        for(j = 0; j < stateDim; j++) {
            double      x_ij = gsl_vector_get(&c_i.vector, j);
            gsl_vector_set(&c_i.vector, j, x_i + s_i * x_ij);
        }
    }

    if( state->lbound!=NULL && state->ubound!=NULL) {
        /*
         * If constraints exist enforce them here. In this case,
         * function values are calculated within constratin_simplex
         * and need not be calculated here.
         */
        status = constrain_simplex(state->x1, state->f1, state->lbound,
                                   state->ubound, func);
        if(status != GSL_SUCCESS)
            GSL_ERROR("constrain_simplex failed", GSL_EFAILED);
    }else{
        /*
         * No constraints exist, so we need to calculate function
         * values.
         */
        for(i = 0; i < stateDim; i++) {
            gsl_vector_view r_i = gsl_matrix_row(&m.matrix, i);
            val = GSL_MULTIMIN_FN_EVAL(minimizer->f, &r_i.vector);
            if(gsl_finite(val)) {
                if(val < state->bestEver) {
                    state->bestEver = val;
                }
            } else {
                GSL_ERROR("non-finite function value encountered",
                          GSL_EBADFUNC);
            }
            gsl_vector_set(state->f1, i + 1, val);
        }
    }

    compute_center(state, state->center);
    (void) compute_size(state, state->center);

    state->bestEver = gsl_vector_min(state->f1);

    sasimplex_sanityCheck(state, __FILE__, __LINE__, __func__);
    return GSL_SUCCESS;
}

/*
 * Do multiple iterations with given temperature. Iterations stop when
 * simplex size declines to tol or when the maximim number, nItr, of
 * iterations is reached. The final simplex size is returned in *size.
 * The function returns the value returned by the final iteration of
 * gsl_multimin_fminimizer_iterate.
 */
int sasimplex_n_iterations(gsl_multimin_fminimizer * minimizer,
                           double *size,
                           double tol_fval,
                           double tol_size,
                           int nItr, double tmptr, int verbose) {
    int         itr = 0, status;

    DPRINTF(("%s:%d:%s: tmptr=%lf\n", __FILE__,__LINE__,__func__,tmptr));

    sasimplex_set_temp(minimizer, tmptr);
    do {
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if(status) {
            printf("%s:%d:%s: %s returned %d: %s\n",
                   __FILE__, __LINE__, __func__,
                   "gsl_multimin_fminimizer_iterate",
                   status, gsl_strerror(status));
            break;
        }

        status = sasimplex_converged(minimizer, size, tol_fval, tol_size);
        switch (status) {
        case GSL_SUCCESS:
            break;
        case GSL_ETOLX:
            /* flat objective function */
            status = GSL_CONTINUE;
            break;
        case GSL_ETOLF:
            /* stuck */
            status = expand_around_best(minimizer->state, minimizer->f);
            if(status != GSL_SUCCESS)
                GSL_ERROR("expand_around_best failed", GSL_EFAILED);
            DPRINTF(("%s:%d:%s: expand_around_best\n",__FILE__,__LINE__,__func__));
            status = GSL_CONTINUE;
            break;
        case GSL_CONTINUE:
            break;
        default:
            fprintf(stderr,
                    "%s:%d:%s: illegal status: %d\n",
                    __FILE__, __LINE__, __func__, status);
            GSL_ERROR("ILLEGAL STATUS", GSL_EFAILED);
        }

        if(verbose) {
            printf("# itr=%d hsiz=%.3f vsiz=%.4f status=%s",
                   itr, *size, sasimplex_vertical_scale(minimizer),
                   statusLbl(status));
            fputs(" x=", stdout);
            vector_print(minimizer->x, stdout);
            putchar('\n');
            fminimizer_print(minimizer);
        }
        ++itr;
    } while(status == GSL_CONTINUE && itr < nItr);

    return status;
}
