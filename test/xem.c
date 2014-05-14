/**
 * @file xem.c
 * @author Alan R. Rogers
 * @brief Test em.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "em.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>


#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         ok = 1, verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xem [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xem [-v]\n");
    }

    unsigned char x, y;

    x = 0;
    y = 3;
    trBits(&x, &y);
    assert(x == 1u && y == 1u);

    x = 3;
    y = 0;
    trBits(&x, &y);
    assert(x == 2u && y == 2u);

    x = 0;
    y = 1;
    trBits(&x, &y);
    assert(x == 0u && y == 1u);

    x = 0;
    y = 0;
    trBits(&x, &y);
    assert(x == 0u && y == 0u);

    x = 2;
    y = 3;
    trBits(&x, &y);
    assert(x == 3 && y == 1);

    x = 2;
    y = 7;
    trBits(&x, &y);
    assert(x == 3u && y == 1u);

    unitTstResult("trBits", "OK");

    x = 0;
    y = 0;
    assert(gamete(x, y) == 0u);

    x = 0;
    y = 1;
    assert(gamete(x, y) == 1u);

    x = 1;
    y = 0;
    assert(gamete(x, y) == 2u);

    x = 1;
    y = 1;
    assert(gamete(x, y) == 3u);

    unitTstResult("gamete", "OK");

    /* EXPERIMENT 1 */

    DsqData     dd = {
        .tol = (double) sqrtl((long double) DBL_EPSILON),
        .px = 0.5,
        .py = 0.5,
        .nGtype = 60,
        .nUnphased = 4,
        .nGam = {29, 29, 29, 29}
    };
    DsqData_reset(&dd);

    if(verbose)
        printf("verbose\n");

    double      lo = loD(dd.px, dd.py, dd.nGam);
    double      hi = hiD(dd.px, dd.py, dd.nGam);

    /* Cause gsl programs to return error codes rather than aborting. */
    gsl_set_error_handler_off();

    if(verbose) {
        printf("Experiment 1: px=%lg py=%lg D in [%lg, %lg]\n",
               dd.px, dd.py, lo, hi);
        printf("  nGam=[%u, %u, %u, %u] nUnphased=%u\n",
               dd.nGam[0], dd.nGam[1], dd.nGam[2], dd.nGam[3], dd.nUnphased);
    }

    double      D;
    int         status;

    status = minimize1D(&D, &dd);
    if(verbose) 
        printf("minimize1D returned status=%d, D=%lg\n", status, D);
    assert(D >= lo && D <= hi);

    /* EXPERIMENT 2 */

    dd.px = 1.0 / 120;
    dd.py = 1.0 / 120;
    dd.nUnphased = 1;
    memset(dd.nGam, 0, 4 * sizeof(dd.nGam[0]));
    dd.nGam[0] = 118;
    DsqData_reset(&dd);

    lo = loD(dd.px, dd.py, dd.nGam);
    hi = hiD(dd.px, dd.py, dd.nGam);

    if(verbose) {
        putchar('\n');
        printf("Experiment 2: px=%lg py=%lg D in [%lg, %lg]\n",
               dd.px, dd.py, lo, hi);
        printf("  nGam=[%u, %u, %u, %u] nUnphased=%u\n",
               dd.nGam[0], dd.nGam[1], dd.nGam[2], dd.nGam[3], dd.nUnphased);
    }

    status = minimize1D(&D, &dd);
    if(verbose)
        printf("minimize1D returned status=%d, D=%lg\n", status, D);
    assert(D >= lo && D <= hi);

    if(verbose)
        putchar('\n');

    /* EXPERIMENT 3 */

    unsigned    twoN = 2 * dd.nGtype;

    dd.px = 20.0 / twoN;
    dd.py = 30.0 / twoN;
    dd.nUnphased = 2;

    dd.nGam[0] = dd.nGam[1] = dd.nGam[2] = dd.nGam[3] = 1;

    lo = loD(dd.px, dd.py, dd.nGam);
    hi = hiD(dd.px, dd.py, dd.nGam);
    D = lo + 0.1 * (hi - lo);

    twoN -= 2 * dd.nUnphased;
    dd.nGam[0] = (unsigned) round(twoN * ((1 - dd.px) * (1 - dd.py) + D));
    dd.nGam[1] = (unsigned) round(twoN * ((1 - dd.px) * dd.py - D));
    dd.nGam[2] = (unsigned) round(twoN * (dd.px * (1 - dd.py) - D));
    dd.nGam[3] = (unsigned) round(twoN * (dd.px * dd.py + D));

    unsigned    twoNtst = dd.nGam[0] + dd.nGam[1] + dd.nGam[2] + dd.nGam[3];

    if(twoN > twoNtst)
        dd.nGam[0] += twoN - twoNtst;

    while(twoN < twoNtst) {
        if(dd.nGam[0])
            --dd.nGam[0];
        else if(dd.nGam[1])
            --dd.nGam[1];
        else if(dd.nGam[2])
            --dd.nGam[2];
        else if(dd.nGam[3])
            --dd.nGam[3];
        twoNtst = dd.nGam[0] + dd.nGam[1] + dd.nGam[2] + dd.nGam[3];
    }
    assert(twoN == dd.nGam[0] + dd.nGam[1] + dd.nGam[2] + dd.nGam[3]);
    DsqData_reset(&dd);

    if(verbose) {
        putchar('\n');
        printf("Experiment 3: px=%lg py=%lg D in [%lg, %lg]\n",
               dd.px, dd.py, lo, hi);
        printf("  nGam=[%u, %u, %u, %u] nUnphased=%u\n",
               dd.nGam[0], dd.nGam[1], dd.nGam[2], dd.nGam[3], dd.nUnphased);
    }

    status = minimize1D(&D, &dd);
    if(verbose)
        printf("minimize1D returned status=%d, D=%lg\n", status, D);
    assert(D >= lo && D <= hi);

    if(verbose)
        putchar('\n');

    /* EXPERIMENT 4: all genotypes are unphased heterozygotes */

    dd.px = 0.5;
    dd.py = 0.5;
    dd.nUnphased = dd.nGtype;

    memset(dd.nGam, 0, sizeof(dd.nGam));

    lo = loD(dd.px, dd.py, dd.nGam);
    hi = hiD(dd.px, dd.py, dd.nGam);

    DsqData_reset(&dd);

    if(verbose) {
        putchar('\n');
        printf("Experiment 4: px=%lg py=%lg D in [%lg, %lg]\n",
               dd.px, dd.py, lo, hi);
        printf("  nGam=[%u, %u, %u, %u] nUnphased=%u\n",
               dd.nGam[0], dd.nGam[1], dd.nGam[2], dd.nGam[3], dd.nUnphased);
    }

    status = minimize1D(&D, &dd);
    if(verbose)
        printf("minimize1D returned status=%d, D=%lg\n", status, D);
    assert(D >= lo && D <= hi);

    if(verbose)
        putchar('\n');

    if(ok)
        unitTstResult("minimize1D", "OK");
    else
        unitTstResult("minimize1D", "FAIL");

    return 0;
}
