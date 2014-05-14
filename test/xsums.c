/**
 * @file xsums.c
 * @author Alan R. Rogers
 * @brief Test sums.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "sums.h"
#include "misc.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#define N 4999
#define TWON (2*(N))

int main(int argc, char **argv) {
    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xsums [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xsums [-v]\n");
    }

    double      t_dotprod, t_dotprod_slow;
    double      t0, t1, t2, t3;
    double      u, v;
    clock_t     start, finish;
    unsigned long i, j, nreps;
    time_t      currtime;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    unsigned long lx[TWON], ly[TWON], ulval;
    unsigned int ix[TWON], iy[TWON], uval;
    unsigned char cx[TWON], cy[TWON];
    unsigned char bcx[N], bcy[N];
    double      x[TWON], y[TWON];

    currtime = time(NULL);
    gsl_rng_set(rng, (unsigned) currtime);

    for(i = 0; i < TWON; ++i) {
        x[i] = gsl_rng_uniform(rng);
        y[i] = gsl_rng_uniform(rng);
        cx[i] = ix[i] = lx[i] = gsl_rng_uniform_int(rng, 2);
        cy[i] = iy[i] = ly[i] = gsl_rng_uniform_int(rng, 2);
    }

    /*
     * bcx and bcy represent each pair of values as a single
     * genotypic value.
     */
    for(i = 0; i < N; ++i) {
        j = 2 * i;
        bcx[i] = cx[j] + 2 * cx[j + 1];
        bcy[i] = cy[j] + 2 * cy[j + 1];
    }

    nreps = 100000;

    assert(sum_long_slow(lx, TWON) == sum_long(lx, TWON));
    unitTstResult("sum_long", "OK");

    u = dotprod_slow(x, y, TWON);
    v = dotprod(x, y, TWON);
    assert(fabs(u - v) <= u * TWON * DBL_EPSILON);
    unitTstResult("dotprod", "OK");

    /* sum_long_slow */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        ulval = sum_long_slow(lx, TWON);
        trickOptimizer();
    }
    finish = clock();
    t0 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose)
        printf("%-16s: %lg sec; val=%lu\n", "sum_long_slow", t0, ulval);

    /* sum_long */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        ulval = sum_long(lx, TWON);
        trickOptimizer();
    }
    finish = clock();
    t1 = ((double) (finish - start)) / CLOCKS_PER_SEC;

    if(verbose) {
        printf("%-16s: %lg sec; val=%lu\n", "sum_long", t1, ulval);
        printf("sum_long speedup: %g%%\n", 100 * (t0 - t1) / t0);
    }

    /* sum_int */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = sum_int(ix, TWON);
        trickOptimizer();
    }
    finish = clock();
    t1 = ((double) (finish - start)) / CLOCKS_PER_SEC;

    if(verbose){
        printf("%-16s: %lg sec; val=%u\n", "sum_int", t1, uval);
        printf("sum_int speedup: %g%%\n", 100 * (t0 - t1) / t0);
    }

    /* sum_char */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = sum_char(cx, TWON);
        trickOptimizer();
    }
    finish = clock();
    t1 = ((double) (finish - start)) / CLOCKS_PER_SEC;

    if(verbose){
        printf("%-16s: %lg sec; val=%u\n", "sum_char", t1, uval);
        printf("sum_char speedup: %g%%\n", 100 * (t0 - t1) / t0);
    }

    /* dotprod_slow */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        dotprod_slow(x, y, TWON);
        trickOptimizer();
    }
    finish = clock();
    t_dotprod_slow = ((double) (finish - start)) / CLOCKS_PER_SEC;

    /* dotprod */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        dotprod(x, y, TWON);
        trickOptimizer();
    }
    finish = clock();
    t_dotprod = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose) {
        printf("dotprod speedup: %lg%%\n",
               100 * (t_dotprod_slow - t_dotprod) / t_dotprod_slow);
    }

    /* dotprod_int_slow */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = dotprod_int_slow(ix, iy, TWON);
        trickOptimizer();
    }
    finish = clock();
    t0 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose)
        printf("%-16s: %lg sec; val=%u\n", "dotprod_int_slow", t0, uval);

    /* dotprod_int */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = dotprod_int(ix, iy, TWON);
        trickOptimizer();
    }
    finish = clock();
    t1 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose)
        printf("%-16s: %lg sec; val=%u\n", "dotprod_int", t1, uval);

    /* sum_and_int */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = sum_and_int(ix, iy, TWON);
        trickOptimizer();
    }
    finish = clock();
    t2 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose) {
        printf("%-16s: %lg sec; val=%u\n", "sum_and_int", t2, uval);

        printf("dotprod_int speedup: %lg%%\n", 100 * (t0 - t1) / t0);
        printf("sum_and_int speedup: %lg%%\n", 100 * (t0 - t2) / t0);
    }

    assert(dotprod_int_slow(ix, iy, TWON) == dotprod_int(ix, iy, TWON));
    unitTstResult("dotprod_int", "OK");

    assert(sum_and_int(ix, iy, TWON) == dotprod_int(ix, iy, TWON));
    unitTstResult("sum_and_int", "OK");

    /* dotprod_char */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = dotprod_char(cx, cy, TWON);
        trickOptimizer();
    }
    finish = clock();
    t1 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose)
        printf("%-16s: %lg sec; val=%u\n", "dotprod_char", t1, uval);

    /* sum_and_char */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = sum_and_char(cx, cy, TWON);
        trickOptimizer();
    }
    finish = clock();
    t2 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose)
        printf("%-16s: %lg sec; val=%u\n", "sum_and_char", t2, uval);

    /* dotprodDiploid */
    start = clock();
    for(i = 0; i < nreps; ++i) {
        uval = dotprodDiploid(bcx, bcy, N);
        trickOptimizer();
    }
    finish = clock();
    t3 = ((double) (finish - start)) / CLOCKS_PER_SEC;
    if(verbose){
        printf("%-16s: %lg sec; val=%u\n", "dotprodDiploid", t3, uval);

        printf("dotprod_char speedup: %lg%%\n", 100 * (t0 - t1) / t0);
        printf("sum_and_char speedup: %lg%%\n", 100 * (t0 - t2) / t0);
        printf("dotprodDiploid speedup: %lg%%\n", 100 * (t0 - t3) / t0);
    }

    assert(dotprod_char(cx, cy, TWON) == dotprod_int(ix, iy, TWON));
    unitTstResult("dotprod_char", "OK");
    
    assert(sum_and_int(ix, iy, TWON) == sum_and_char(cx, cy, TWON));
    unitTstResult("sum_and_char", "OK");

    assert(sum_and_char(cx, cy, TWON) == dotprodDiploid(bcx, bcy, N));
    unitTstResult("dotprodDiploid", "OK");

    assert(sum_char(cx, TWON) == sumDiploid(bcx, N));
    unitTstResult("sumDiploid", "OK");

    gsl_rng_free(rng);

    return (0);
}
