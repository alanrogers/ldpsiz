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

const mpfr_rnd_t rnd = MPFR_RNDN;  // round to nearest
const mpfr_prec_t precision = 128; // bits of precision

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
    mpfr_t *rvec; // row eigenvectors
    mpfr_t *cvec; // column eigenvectors
    mpfr_t *beta; // beta[i-2] = i*(i-1)/2
    mpfr_t *lambda; // eigenvalues
    mpfr_t *x;  // work vector
    unsigned *offset;  // offset[j] = j*(j+1)/2
};

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim]);

/// Allocate and initialize a new vector of mpfr_t values.
/// If x==NULL on input, the vector is allocated but not initialized.
MpfrVec *MpfrVec_new(unsigned dim, long double x[dim]) {
    MpfrVec *new = malloc(sizeof(new[0]));
    CHECKMEM(new);

    new->dim = dim;
    new->x = malloc(dim * sizeof(new->x[0]));
    unsigned i;
    for(i=0; i < dim; ++i)
        mpfr_init2(new->x[i], precision);

    if(x != NULL)
        for(i=0; i<dim; ++i)
            mpfr_set_ld(new->x[i], x[i], rnd);

    return new;
}

/// Deallocate MpfrVec
void MpfrVec_free(MpfrVec *self) {
    unsigned i;
    for(i=0; i < self->dim; ++i)
        mpfr_clear(self->x[i]);
    free(self->x);
    free(self);
}

/// Copy contents of MpfrVec into array x.
void MpfrVec_get(MpfrVec *self, unsigned dim, long double x[dim]) {
    assert(dim == self->dim);
    unsigned i;

    for(i=0; i < dim; ++i)
        x[i] = mpfr_get_ld(self->x[i], rnd);
}

MatCoalSpec *MatCoalSpec_new(unsigned nSamples) {
    long i, j, ii, jj;
    mpfr_t x, y, z;
    MatCoalSpec *self = malloc( sizeof *self );
    CHECKMEM(self);

    self->nSamples = nSamples;
    self->dim = nSamples-1;
    self->nPairs = (self->dim * (self->dim+1))/2;

    mpfr_init2(x, precision);
    mpfr_init2(y, precision);
    mpfr_init2(z, precision);

    self->offset = malloc(self->dim * sizeof self->offset[0]);
    CHECKMEM(self->offset);

    self->rvec = malloc(self->nPairs * sizeof self->rvec[0]);
    CHECKMEM(self->rvec);

    self->cvec = malloc(self->nPairs * sizeof self->cvec[0]);
    CHECKMEM(self->cvec);

    for(i=0; i < self->nPairs; ++i) {
        mpfr_init2(self->rvec[i], precision);
        mpfr_init2(self->cvec[i], precision);
    }

    for(i=0; i < self->dim; ++i)
        self->offset[i] = (i*(i + 1))/2;

    for(j=2; j <= nSamples; ++j) {
        jj = j-2;
        mpfr_set_si(self->cvec[jj + self->offset[jj]], 1L, rnd);
        mpfr_set_si(self->rvec[jj + self->offset[jj]], 1L, rnd);
        for(i = j-1; i > 1; --i) {
            ii = i-2;
            mpfr_set_si(x, i*(i+1L), rnd);
            mpfr_set_si(y, i*(i-1L) - j*(j-1L), rnd);
            mpfr_div(z, x, y, rnd);
            // now z = i*(i+1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->cvec[ii + self->offset[jj]],
                     self->cvec[ii+1 + self->offset[jj]],
                     z, rnd);
        }
        for(i=j+1; i <= nSamples; ++i) {
            ii = i-2;
            mpfr_set_si(x, i*(i-1L), rnd);
            mpfr_set_si(y, i*(i-1L) - j*(j-1L), rnd);
            mpfr_div(z, x, y, rnd);
            // z = i*(i-1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->rvec[jj+self->offset[ii]],
                     self->rvec[jj+self->offset[ii-1]], z, rnd);
        }
    }

    self->beta = malloc(self->dim * sizeof self->beta[0]);
    CHECKMEM(self->beta);

    self->lambda = malloc(self->dim * sizeof self->lambda[0]);
    CHECKMEM(self->lambda);

    self->x = malloc(self->dim * sizeof self->x[0]);
    CHECKMEM(self->x);

    // eigenvalues are beta[i]*t/N
    for(i=0; i < self->dim; ++i) {
        j = i+2;
        mpfr_init2(self->beta[i], precision);
        mpfr_init2(self->lambda[i], precision);
        mpfr_init2(self->x[i], precision);
        mpfr_set_si(self->beta[i], (j*(j-1L))/2L, rnd);
    }
    mpfr_clears(x, y, z, (mpfr_ptr) 0);
    return self;
}

void MatCoalSpec_project(MatCoalSpec *self, MpfrVec *x, long double v) {
    assert(x);
    if(x->dim != self->dim)
        eprintf("%s:%s:%d: dimensions don't match. x->dim=%ld but self->dim=%u\n",
                __FILE__,__func__,__LINE__, x->dim, self->dim);
    int i;
    mpfr_t work;
    mpfr_init2(work, precision);
    
    // Left-multiply state vector, x, by matrix of row eigenvectors.
    UTmatXvec(self->dim, self->x, self->rvec, self->offset, x->x);

    // Left-multiply resulting vector by diag(exp(-beta[i]*v))
    for(i=0; i < self->dim; ++i) {
        // work = exp(-beta[i]*v)
        mpfr_set_d(work, v, rnd);
        mpfr_mul(work, work, self->beta[i], rnd); 
        mpfr_neg(work, work, rnd); 
        mpfr_exp(work, work, rnd); 
        mpfr_mul(self->x[i], self->x[i], work, rnd);
    }

    // Left-multiply by matrix of column eigenvectors.
    UTmatXvec(self->dim, x->x, self->cvec, self->offset, self->x);

    mpfr_clear(work);
}

/// Form matrix product y = A*x, where y is a vector, A an
/// upper-triangular matrix, and x is a vector.  x and y have
/// dimension dim, and A has dimension dim X dim. offset[j] is the
/// offset of the beginning of the j'th column within matrix A.
/// A should be laid out so that element (i,j), i.e. row i and column
/// j, is at position A[i + offset[j]].
void UTmatXvec(unsigned dim, mpfr_t *y, mpfr_t *A, unsigned *offset, mpfr_t *x) {
    unsigned i, j;
    mpfr_t tmp;
    mpfr_init2(tmp, precision);

    for(i=0; i<dim; ++i) {
        mpfr_set_d(y[i], 0, rnd);
        for(j=i; j<dim; ++j) {
            mpfr_mul(tmp, A[i + offset[j]], x[j], rnd);
            mpfr_add(y[i], y[i], tmp, rnd);
        }
    }

    mpfr_clear(tmp);
}

void MatCoalSpec_print(MatCoalSpec *self) {
    printf("nSamples=%u nPairs=%d precision=%ld\n",
           self->nSamples, self->nPairs, precision);

    long i;

    printf("beta:");
    for(i=0; i < self->dim; ++i)
        mpfr_printf(" %RNf", self->beta[i]);
    putchar('\n');
    
    printf("Row eigenvectors:\n");
    prUTMat(self->dim, self->rvec, 17, self->offset);

    printf("Column eigenvectors:\n");
    prUTMat(self->dim, self->cvec, 17, self->offset);
}

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim]) {
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

    unsigned i, nSamples = 5;
    MatCoalSpec *mcs = MatCoalSpec_new(nSamples);

    if(verbose)
        MatCoalSpec_print(mcs);

    unsigned  dim = nSamples-1;
    long double x[dim], y[dim];
    memset(x, 0, sizeof(x));
    x[dim-1] = 1.0L;

    MpfrVec *vec = MpfrVec_new(dim, x);

    long double v = 0.0L;
    MatCoalSpec_project(mcs, vec, (long double) 0.0L);
    MpfrVec_get(vec, dim, y);
    printf("v=%Lf\n", v);
    for(i=0; i<dim; ++i)
        printf("%d %10.6Lf -> %10.6Lf\n", i+2, x[i], y[i]);

    v = 1.0L;
    MatCoalSpec_project(mcs, vec, v);
    MpfrVec_get(vec, dim, y);
    printf("v=%Lf\n", v);
    for(i=0; i<dim; ++i)
        printf("%d %10.6Lf -> %10.6Lf\n", i+2, x[i], y[i]);

    v = 100000.0L;
    MatCoalSpec_project(mcs, vec, v);
    MpfrVec_get(vec, dim, y);
    printf("v=%Lf\n", v);
    for(i=0; i<dim; ++i)
        printf("%d %10.6Lf -> %10.6Lf\n", i+2, x[i], y[i]);

    return 0;
}
#endif
