/**
 * @file sums.c
 * @author Alan R. Rogers
 * @brief Optimized functions for sums and dot products.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <assert.h>
#include "misc.h"
#include "sums.h"

/*
 * Form dot product (i.e. inner product) of vectors, x and y, of doubles. 
 *
 * Alan Rogers 2012-4-10.
 */
double dotprod(double *x, double *y, unsigned n) {
    register double rval;
    register unsigned i, m;

    myassert(n >= 0);

    rval = 0.0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i] * y[i];

    for(i = m; i < n; i += 5) {
        rval += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return rval;
}

/*
 * Form dot product (i.e. inner product) of unsigned vectors x and y. 
 *
 * Alan Rogers 2012-4-10.
 */
unsigned dotprod_int(unsigned *x, unsigned *y, unsigned n) {
    register unsigned rval, i, m;

    myassert(n >= 0);

    rval = 0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i] * y[i];

    for(i = m; i < n; i += 5)
        rval += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];

    return rval;
}

/*
 * Form dot product (i.e. inner product) of char vectors x and y. 
 *
 * Alan Rogers 2012-4-10.
 */
unsigned dotprod_char(unsigned char *x, unsigned char *y, unsigned n) {
    register unsigned rval, i, m;

    myassert(n >= 0);

    rval = 0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i] * y[i];

    for(i = m; i < n; i += 5)
        rval += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];

    return rval;
}

/*
 * Form dot product (i.e. inner product) of vectors x and y, in which
 * each char represents a 2-bit binary number.
 *
 * At each locus, there are 2 alleles (0 and 1) and 4 ordered
 * genotypes: 00, 01, 10, and 11. These are binary integers with
 * decimal values 0, 1, 2, and 3. The entries of x and y are treated
 * as numbers in this format. Before calling this function, the user
 * should ensure that each value is in the range [0,3].
 *
 * @param[in] x,y Arrays of ordered diploid genotypes, each represented as a
 * binary number in the range [0,3].
 *
 * @param[in] n The number of diploid genotypes in each array.
 *
 * @returns The inner product of the two arrays. Here, multiplication
 * is defined to give the number of '1' bits in the binary "and" of
 * the two values. This is equivalent to the cross product of x and y
 * across all 2n bits.
 *
 * Alan Rogers 2012-12-21.
 */
unsigned dotprodDiploid(unsigned char *x, unsigned char *y, unsigned n) {
    register unsigned rval, i, k;

    /* multiplication matrix */
    static const unsigned char M[4][4] = { {0, 0, 0, 0},
    {0, 1, 0, 1},
    {0, 0, 1, 1},
    {0, 1, 1, 2}
    };

    myassert(n >= 0);

    rval = 0;
    k = n % 5;
    for(i = 0; i < k; ++i) {
        assert(x[i] < 4);
        assert(y[i] < 4);
        rval += M[x[i]][y[i]];
    }

    for(i = k; i < n; i += 5)
        rval += M[x[i]][y[i]]
            + M[x[i + 1]][y[i + 1]]
            + M[x[i + 2]][y[i + 2]]
            + M[x[i + 3]][y[i + 3]]
            + M[x[i + 4]][y[i + 4]];

    return rval;
}

unsigned sumDiploid(const unsigned char *x, unsigned n) {
    register unsigned rval, i, k;

    /*
     * This function assumes that UNPHASED_HETEROZYGOTE==4.
     */
    assert(UNPHASED_HETEROZYGOTE == 4);

    /* Number of copies of "1" allele in each genotype */
    static const unsigned char S[5] = { 0, 1, 1, 2, 1 };

    myassert(n >= 0);

    rval = 0;
    k = n % 5;
    for(i = 0; i < k; ++i) {
        assert(x[i] < 5);
        rval += S[x[i]];
    }

    for(i = k; i < n; i += 5)
        rval += S[x[i]]
            + S[x[i + 1]]
            + S[x[i + 2]]
            + S[x[i + 3]]
            + S[x[i + 4]];

    return rval;
}

/*
 * Sum the values of x[i] & y[i]. If x and y are binary (0/1) vectors,
 * this yields the number of pairs of form (1,1) and is equivalent to
 * the dot product.
 *
 * Alan Rogers 2012-4-11.
 */
unsigned sum_and_int(unsigned *x, unsigned *y, unsigned n) {
    register unsigned rval, i, m;

    myassert(n >= 0);

    rval = 0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i] & y[i];

    for(i = m; i < n; i += 5) {
        rval += (x[i] & y[i]) + (x[i + 1] & y[i + 1]) + (x[i + 2] & y[i + 2])
            + (x[i + 3] & y[i + 3]) + (x[i + 4] & y[i + 4]);
    }

    return rval;
}

unsigned sum_and_char(unsigned char *x, unsigned char *y, unsigned n) {
    register unsigned rval, i, m;

    myassert(n >= 0);

    rval = 0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i] & y[i];

    for(i = m; i < n; i += 5) {
        rval += (x[i] & y[i]) + (x[i + 1] & y[i + 1]) + (x[i + 2] & y[i + 2])
            + (x[i + 3] & y[i + 3]) + (x[i + 4] & y[i + 4]);
    }

    return rval;
}

unsigned long sum_long(unsigned long *x, unsigned n) {
    register unsigned long rval;
    register unsigned i, m;

    myassert(n >= 0);

    rval = 0L;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i];

    for(i = m; i < n; i += 5)
        rval += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];

    return rval;
}

double sum_double(double *x, unsigned n) {
    register double rval;
    register unsigned i, m;

    myassert(n >= 0);

    rval = 0.0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i];

    for(i = m; i < n; i += 5)
        rval += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];

    return rval;
}

unsigned sum_int(unsigned *x, unsigned n) {
    register unsigned rval;
    register unsigned i, m;

    myassert(n >= 0);

    rval = 0.0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i];

    for(i = m; i < n; i += 5)
        rval += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];

    return rval;
}

unsigned sum_char(const unsigned char *x, unsigned n) {
    register unsigned rval;
    register unsigned i, m;

    myassert(n >= 0);

    rval = 0.0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i];

    for(i = m; i < n; i += 5)
        rval += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];

    return rval;
}

double dotprod_slow(double *x, double *y, unsigned n) {
    register double rval;
    register unsigned i;

    myassert(n >= 0);

    rval = 0.0;
    for(i = 0; i < n; ++i)
        rval += x[i] * y[i];

    return rval;
}

unsigned long sum_long_slow(unsigned long *x, unsigned n) {
    register unsigned long rval;
    register unsigned i;

    myassert(n >= 0);

    rval = 0.0;
    for(i = 0; i < n; ++i)
        rval += x[i];

    return rval;
}

unsigned dotprod_int_slow(unsigned *x, unsigned *y, unsigned n) {
    register unsigned rval, i;

    myassert(n >= 0);

    rval = 0;
    for(i = 0; i < n; ++i)
        rval += x[i] * y[i];

    return rval;
}

unsigned dotprod_char_slow(unsigned char *x, unsigned char *y, unsigned n) {
    register unsigned rval, i;

    myassert(n >= 0);

    rval = 0;
    for(i = 0; i < n; ++i)
        rval += x[i] * y[i];
    return rval;
}

void trickOptimizer(void) {
}
