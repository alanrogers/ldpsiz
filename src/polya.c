/**
   @file polya.c
   @author Alan R. Rogers
   @brief Miscellaneous functions.
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "polya.h"
#include "typedefs.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static void polyaVector(double  *p, int n, int k);

/// Constructor
/// @param n The number of genes in the sample
Polya *Polya_new(int n) {
    Polya *polya = malloc(sizeof(Polya));
    checkmem(polya, __FILE__, __LINE__);

    polya->n = n;
    polya->mat = malloc((n-1)*n*sizeof(polya->mat[0]));
    checkmem(polya->mat, __FILE__, __LINE__);

    int k;
    for(k=2; k <= n; ++k)
        polyaVector(polya->mat + (k-2)*n, n, k);
    return polya;
}

void Polya_free(Polya *polya) {
    assert(polya);
    free(polya->mat);
    free(polya);
}

/// Print Polya distribution
void Polya_print(const Polya *polya, FILE *ofp) {
    int i, k;

    printf("Pr[i,k]=Prob i mutants in sample of %d,"
           " given mutation occurred in interval k\n",  polya->n);
    printf("%4s", "i\\k");
    for(k=2; k <= polya->n; ++k) 
        printf(" %8d", k);
    putchar('\n');
    for(i=1; i < polya->n; ++i) {
        printf("%4d", i);
        for(k=2; k <= polya->n; ++k) 
            printf(" %8.4lf", Polya_prob(polya, i, k));
        putchar('\n');
    }
}

/**
   @param p vector of length at least n-k+1.  On return, p[i] =
   Pr[observing i+1 mutants in a sample of size n, given that the
   mutation occurred at a time when the gene tree had k lineages].

   @param n sample size

   @param k the number of lineages in the gene tree when the mutation
   occurred. 
**/
static void polyaVector(double  *p, int n, int k) {

    int i, nn;

	assert(k >= 2);
	assert(k <= n);
    if(n < 2)
        eprintf("%s:%s:%d: illegal argument n=%d\n",
                __FILE__, __func__, __LINE__, n);

	// polya distribution for nn=k
	p[0] = 1;
	p[1] = 0;   // writes outside array if n<2.

	// nn represents sample size
	for(nn=k+1; nn <= n; ++nn) {
	    double a=0;
	    double b=p[0];

	    for(i=1; ; ++i) {
            double f = ((double)(i-1))/(nn-1);
            double g = 1.0 - ((double)i)/(nn-1);
            p[i-1] = a*f + b*g;
            if( i == nn-k+1 ) {
                if(i < n-1)
                    // max i = n-2 if k == 2,3
                    //       = n-k+1 if k > 4,5,...
                    p[i] = 0;
                break;
            }
            a = b;
            b = p[i];
	    }

	    // For k=2,3 we have set p[0]..p[n-2].  For k=4,5,... we have
	    // set p[0]...p[n-k+1].  In the latter case, fill
	    // p[n-k+2]...p[n-2] with zeroes.  This amounts to an
	    // additional k-3 zeroes. 
	    if(k>3)
            memset(p+n-k+2, 0, (k-3) * sizeof(double));
	}
}




