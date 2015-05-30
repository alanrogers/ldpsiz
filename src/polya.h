/**
   @file polya.h
   @brief Header for polya.c
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#ifndef MATCOAL_POLYA_H
#define MATCOAL_POLYA_H

#include "typedefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct Polya {
    // Number of genes in sample
    int n;

    // mat represents a matrix of n-1 rows and n columns. But this
    // matrix is represented here as a 1-dimensional array. 
    //
    // Row k contains the polya distribution for mutations that
    // occur in a coalescent interval containing k+1 lines of descent.
    //
    // Column i is the probability that such a mutation has i+i descendants
    // in a sample of size n. This probability is stored as mat[k*n + i]. 
    double *mat; 
};

Polya *Polya_new(int n);
void Polya_free(Polya *polya);
void Polya_print(const Polya *polya, FILE *ofp);
static inline double Polya_prob(const Polya *polya, int i, int k);

/// Probability of i mutants in sample given that mutation occurred
/// in a coalescent interval containing k lineages.
static inline double Polya_prob(const Polya *polya, int i, int k) {
    assert(polya);
#ifndef NDEBUG
    if(k < 2 || k > polya->n) {
        fprintf(stderr,"%s:%d: invalid value: k=%d; must satisfy 2 < k <= %d\n",
                __FILE__, __LINE__, k, polya->n);
        exit(1);
    }
#endif
	    
    if( i < 1 || i > polya->n - k + 1)
        return 0;

    return polya->mat[(k-2)*polya->n + i-1];
}

#endif
