/**
   @file xpolya.c
   @brief Test polya.c
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "polya.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

int main(int argc, char **argv) {

	int i, n=7, ok=1, verbose=0;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose = 1;
        else 
            eprintf("usage: xpolya [-v]\n");
    }

    double calc, truth, err, errTol = 0.00001;
	Polya *polya = Polya_new(n);
	checkmem(polya, __FILE__, __LINE__);

    if(verbose)
        Polya_print(polya, stdout);

    if(Polya_prob(polya, 1,n) != 1) {
        ok = 0;
        if(verbose)
            printf("Polya_prob(polya, 1,n) = %lg; should = 1; FAIL\n",
                   Polya_prob(polya, 1,n));
    }
    if(Polya_prob(polya, 2,n) != 0) {
        ok = 0;
        if(verbose)
            printf("Polya_prob(polya, 2,n) = %lg; should = 0; FAIL\n",
                   Polya_prob(polya, 2,n));
    }
	for(i=1; i < n; ++i) {
        calc = Polya_prob(polya,i,2);
        truth = 1.0/(n-1);
        err = fabs(calc-truth);
        if(verbose || err>errTol)
            printf("Polya(%d,2)=%lg truth=%lg err=%lg",
                   i, calc, truth, err);
        if(err > errTol) {
            ok=0;
            printf(" FAIL\n");
        }else if(verbose)
            putchar('\n');
	}
    i=3;
    calc = Polya_prob(polya,3,4);
    truth = 0.15;
    err = fabs(calc-truth);
    if(verbose || err>errTol)
        printf("Polya(%d,4)=%lg truth=%lg err=%lg",
               i, calc, truth, err);
    if(err > errTol) {
        ok=0;
        printf(" FAIL\n");
    }else if(verbose)
        putchar('\n');
    
	Polya_free(polya);

	unitTstResult("Polya", ok ? "OK" : "FAIL");
}
