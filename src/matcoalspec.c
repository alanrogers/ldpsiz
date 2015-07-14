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
    unsigned dim; // dimension of matrices and vectors
    unsigned nPairs;
    mpfr_t *rvec; // row eigenvectors
    mpfr_t *cvec; // column eigenvectors

    // offset[j] = j*(2K-j+1)/2
    unsigned *offset;
};

MatCoalSpec *MatCoalSpec_new(unsigned dim) {
    MatCoalSpec *new = malloc( sizeof *new );
    CHECKMEM(new);

    new->dim = dim;
    new->nPairs = (dim*(dim-1))/2;

    new->offset = malloc(dim * sizeof new->offset[0]);
    CHECKMEM(new->offset);

    new->rvec = malloc(new->nPairs * sizeof new->rvec[0]);
    CHECKMEM(new->rvec);

    new->cvec = malloc(new->nPairs * sizeof new->cvec[0]);
    CHECKMEM(new->cvec);

    unsigned i, j, ii, jj;
    for(i=0; i < dim; ++i)
        new->offset[i] = (i*(2*dim - i + 1))/2;

    XXXXXXXXXXXXXXXXX stopped here
}


