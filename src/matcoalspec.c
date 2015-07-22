/**
 * @file matcoalspec.c
 * @brief Spectral decomposition of matrix coalescent
 * @internal
 * Copyright (C) 2015 Alan R. Rogers
 *
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * This file is released under the Internet Systems Consortium
 * License, which can be found in file "LICENSE".
 *
 * Alan R. Rogers, Department of Anthropology, University of Utah,
 * Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
 **/

#include "matcoalspec.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <mpfr.h>
#include <math.h>

/*
   0  1  2  3  4  col j
   0  1  3  6 10  offset = j*(j+1)/2
  00 01 02 03 04  
  xx 11 12 13 14  
  xx xx 22 23 24  
  xx xx xx 33 34  
  xx xx xx xx 44  

  offset[j] = j*(j+1)/2
  x[i][j] is array[offset[j] + i]
 */
struct MatCoalSpec {
    unsigned nSamples;
    unsigned dim; // dimension of matrices and vectors
    unsigned nPairs;
    mpfr_prec_t precision; // binary(?) digits of precision
    mpfr_rnd_t rnd;  // rounding mode
    mpfr_t *rvec; // row eigenvectors
    mpfr_t *cvec; // column eigenvectors
    mpfr_t *beta; // beta[i-2] = i*(i-1)/2
    mpfr_t *lambda; // eigenvalues
    mpfr_t *x;      // x[i] = probability of state i+2
    unsigned *offset;  // offset[j] = j*(j+1)/2
};

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim],
             mpfr_rnd_t rnd);

MatCoalSpec *MatCoalSpec_new(unsigned nSamples, unsigned precision) {
    long i, j, ii, jj;
    mpfr_t x, y, z;
    MatCoalSpec *self = malloc( sizeof *self );
    CHECKMEM(self);

    self->nSamples = nSamples;
    self->dim = nSamples-1;
    self->nPairs = (self->dim * (self->dim+1))/2;
    self->precision = precision;
    self->rnd = MPFR_RNDN; // round to nearest

    mpfr_init2(x, self->precision);
    mpfr_init2(y, self->precision);
    mpfr_init2(z, self->precision);

    self->offset = malloc(self->dim * sizeof self->offset[0]);
    CHECKMEM(self->offset);

    self->rvec = malloc(self->nPairs * sizeof self->rvec[0]);
    CHECKMEM(self->rvec);

    self->cvec = malloc(self->nPairs * sizeof self->cvec[0]);
    CHECKMEM(self->cvec);

    for(i=0; i < self->nPairs; ++i) {
        mpfr_init2(self->rvec[i], self->precision);
        mpfr_init2(self->cvec[i], self->precision);
    }

    for(i=0; i < self->dim; ++i)
        self->offset[i] = (i*(i + 1))/2;

    for(j=2; j <= nSamples; ++j) {
        jj = j-2;
        mpfr_set_si(self->cvec[jj + self->offset[jj]], 1L, self->rnd);
        mpfr_set_si(self->rvec[jj + self->offset[jj]], 1L, self->rnd);
        for(i = j-1; i > 1; --i) {
            ii = i-2;
            mpfr_set_si(x, i*(i+1L), self->rnd);
            mpfr_set_si(y, i*(i-1L) - j*(j-1L), self->rnd);
            mpfr_div(z, x, y, self->rnd);
            // now z = i*(i+1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->cvec[ii + self->offset[jj]],
                     self->cvec[ii+1 + self->offset[jj]],
                     z, self->rnd);
        }
        for(i=j+1; i <= nSamples; ++i) {
            ii = i-2;
            mpfr_set_si(x, i*(i-1L), self->rnd);
            mpfr_set_si(y, i*(i-1L) - j*(j-1L), self->rnd);
            mpfr_div(z, x, y, self->rnd);
            // z = i*(i-1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->rvec[jj+self->offset[ii]],
                     self->rvec[jj+self->offset[ii-1]], z, self->rnd);
        }
    }

    self->beta = malloc(self->dim * sizeof self->beta[0]);
    CHECKMEM(self->beta);

    self->lambda = malloc(self->dim * sizeof self->lambda[0]);
    CHECKMEM(self->beta);

    self->x = malloc(self->dim * sizeof self->x[0]);
    CHECKMEM(self->beta);

    // eigenvalues are beta[i]*t/N
    for(i=0; i < self->dim; ++i) {
        j = i+2;
        mpfr_init2(self->beta[i], self->precision);
        mpfr_init2(self->lambda[i], self->precision);
        mpfr_init2(self->x[i], self->precision);
        mpfr_set_si(self->beta[i], (j*(j-1L))/2L, self->rnd);
    }

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(z);

    return self;
}

int MatCoalSpec_project(MatCoalSpec *self, int dim,
                        long double x[dim], long double v) {
    if(dim != self->dim)
        eprintf("%s:%s:%d: dimensions don't match. dim=%ld but self->dim=%u\n",
                __FILE__,__func__,__LINE__, dim, self->dim);
    int i;
    mpfr_t lambda, v2;
    mpfr_init2(lambda, self->precision);
    mpfr_init2(v2, self->precision);
    mpfr_set_d(v2, v, self->rnd);
    if(x != NULL) {
        // use x vector passed as argument
        for(i=0; i < dim; ++i)
            mpfr_set_ld(self->y[i], x[i], self->rnd);
    }else {
        // use x vector from previous call to MatCoalSpec_project
        for(i=0; i < dim; ++i)
            mpfr_set(self->y[i], self->x[i], self->rnd);
    }

    // Left-multiply state vector, y, by matrix of row eigenvectors.
    UTmatXvec(self->dim, self->x, self->rvec, self->y, self->rnd);

    // Left-multiply resulting vector by diag(exp(-beta[i]*v))
    for(i=0; i < dim; ++i) {
        // lambda[i] = exp(-beta[i]*v)
        mpfr_mul(lambda, self->beta[i], v2); 
        mpfr_neg(lambda, v2); 
        mpfr_exp(lambda, v2); 
        mpfr_mul(self->y[i], self->x[i], lambda, v2); 
    }

    // Left-multiply by matrix of column eigenvectors.
    UTmatXvec(self->dim, self->x, self->rvec, self->y, self->rnd);

    // Store result in x.
    for(i=0; i < dim; ++i)
        x[i] = mpfr_get_d(self->x[i], self->rnd);
}

void MatCoalSpec_print(MatCoalSpec *self) {
    printf("nSamples=%u nPairs=%d precision=%ld\n",
           self->nSamples, self->nPairs, self->precision);

    long i;

    printf("beta:");
    for(i=0; i < self->dim; ++i)
        mpfr_printf(" %RNf", self->beta[i]);
    putchar('\n');
    
    printf("Row eigenvectors:\n");
    prUTMat(self->dim, self->rvec, 17, self->offset, self->rnd);

    printf("Column eigenvectors:\n");
    prUTMat(self->dim, self->cvec, 17, self->offset, self->rnd);
}

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim],
             mpfr_rnd_t rnd) {
    int i, j;
    for(i=0; i<dim; ++i) {
        for(j=0; j < i; ++j) 
            printf(" %*d", prWid, 0);
        for(j=i; j < dim; ++j) {
            putchar(' ');
#if 0
            cWritten = mpfr_out_str(stdout, 10, prWid,
                                    mat[i + offset[j]],
                                    rnd);
#else
            mpfr_printf("%*.6RNf", prWid,mat[i + offset[j]]);
#endif
        }
        putchar('\n');
    }
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xmatcoalspec [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    unsigned nSamples = 5;
    unsigned precision = 128;

    MatCoalSpec *mcs = MatCoalSpec_new(nSamples, precision);

    if(verbose)
        MatCoalSpec_print(mcs);

    return 0;
}
#endif
