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
#include <stdarg.h>
#include <mpfr.h>
#include <math.h>

/*                                           j*(2K-j+1)/2
  00 01 02 03 04  row 0: offset 0                       0
  xx 11 12 13 14  row 1: offset 5                       5
  xx xx 22 23 24  row 2: offset 5+4=9                   9
  xx xx xx 33 34  row 3: offset 5+4+3=12               12
  xx xx xx xx 44  row 4: offset 5+4+3+2=14             14

  offset[j] = j*(2K-j+1)/2
 */
struct MatCoalSpec {
    unsigned nSamples;
    unsigned dim; // dimension of matrices and vectors
    unsigned nPairs;
    mpfr_prec_t precision; // binary(?) digits of precision
    mpfr_rnd_t rnd;  // rounding mode
    mpfr_t *rvec; // row eigenvectors
    mpfr_t *cvec; // column eigenvectors

    // offset[j] = j*(2K-j+1)/2
    unsigned *offset;
};

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim],
             mpfr_rnd_t rnd);

MatCoalSpec *MatCoalSpec_new(unsigned nSamples, unsigned precision) {
    unsigned i, j, ii, jj;
    mpfr_t lambda; // eigenvalue
    mpfr_t x;
    MatCoalSpec *self = malloc( sizeof *self );
    CHECKMEM(self);

    self->nSamples = nSamples;
    self->dim = nSamples-2;
    self->nPairs = (self->dim * (self->dim-1))/2;
    self->precision = precision;
    self->rnd = MPFR_RNDN; // round to nearest

    mpfr_init2(lambda, self->precision);
    mpfr_init2(x, self->precision);

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
        self->offset[i] = (i*(2 * self->dim - i + 1))/2;

    for(j=2; j < nSamples; ++j) {
        unsigned ii, jj;
        long double tmp;
        jj = j-2;
        long double lambda = -j*(j-1.0);
        mpfr_set_d(self->cvec[self->offset[jj]+jj], 1.0, self->rnd);
        for(i = j-1; i > 1; --i) {
            ii = i-2;
            tmp = i*(i+1.0L)/(i*(i-1.0L) + lambda);
            mpfr_set_ld(x, tmp, self->rnd);
            mpfr_mul(self->cvec[self->offset[ii]+jj],
                     self->cvec[self->offset[ii+1]+jj], x, self->rnd);
        }
        for(i=j+1; i <= nSamples; ++i) {
            ii = i-2;
            tmp = i*(i-1.0)/(i*(i-1.0) + lambda);
            mpfr_set_ld(x, tmp, self->rnd);
            rvec[jj][ii] = rvec[jj][ii-1]*tmp;
            mpfr_mul(self->rvec[self->offset[jj]+ii],
                     self->rvec[self->offset[jj]+ii-1], x, self->rnd);
        }
    }

    size_t cWritten;

    printf("Row eigenvectors:\n");
    prUTMat(self->dim, self->rvec, 8, self->offset, self->rnd);

    printf("Column eigenvectors:\n");
    prUTMat(self->dim, self->cvec, 8, self->offset, self->rnd);

    mpfr_clear(lambda);
    mpfr_clear(x);
}

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim],
             mpfr_rnd_t rnd) {
    int i, j;

    for(i=0; i<dim; ++i) {
        for(j=0; j < i; ++j) 
            printf(" %*d", prWid, 0);
        for(j=i; j < i; ++j) {
            putchar(' ');
            cWritten = mpfr_out_str(stdout, 10, prWid,
                                    mat[offset[i]+j],
                                    rnd);
        }
        putchar('\n');
    }
}


