/**
 * @file em.h
 * @author Alan R. Rogers
 * @brief Header for em.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_EM_H
#define LDPSIZ_EM_H

#include <stdio.h>
#include <assert.h>
#include "typedefs.h"
#include "window.h"

/** Is A a homozygote? */
#define HOMOZYGOTE(A) ((A)==0 || (A)==3)

/** True if A and B are both heterozygotes and at least one is unphased. */
#define UNPHASED_PAIR(A,B) (((A)==UNPHASED_HETEROZYGOTE && !HOMOZYGOTE(B)) \
                            || ((B)==UNPHASED_HETEROZYGOTE && !HOMOZYGOTE(A)))

/**
 * If A is an unphased heterozygote, then convert it to a phased
 * one (genotype 1). Otherwise leave it alone.
 */
#define MAKE_PHASED(A) ((A)==UNPHASED_HETEROZYGOTE ? 1 : (A))

/**
 * Data used by EM algorithm in estimating Dsq from partially phased
 * diploid data.
 */
typedef struct DsqData {
    double      tol;        /**< controls convergence */
    double      alpha;      /**< a*(1-a)*b*(1-b) */
    double      beta;       /**< a + b - 2*a*b */
    double      loD;        /**< least feasible D */
    double      hiD;        /**< greatest feasible D */
    double      px;         /**< frequency of allele "1" in sample from A */
    double      py;         /**< frequency of allele "1" in sample from B */
    unsigned    nGtype;     /**< number of genotypes in (diploid) sample */
    unsigned    nUnphased;  /**< number of unphased genotype pairs */
    unsigned    nGam[4];    /**< number of phased gametes of types 0-3 */
} DsqData;

void        DsqData_print(DsqData * dd, const char *file, int line,
                          FILE * fp);
void        DsqData_reset(DsqData * dd);
double      loD(double pA, double pB, unsigned *nGam);
double      hiD(double pA, double pB, unsigned *nGam);
int         minimize1D(double *D, DsqData * dd);

/**
 * Think of x and y as the two rows of a matrix, with the columns of
 * the matrix equal to the bits of x and y. We are only interested in
 * the low-order (rightmost) pair of bits. Transpose that 2X2 matrix,
 * so that x and y are set equal to the bits that had previously been
 * in the columns.
 *
 * In other words, if [i,j] are the right-most two bits of x, and
 * [k,l] are the right-most two bits of y, then the function does the
 * following transformation:
 *
 *  x:  i  j             x:  i  k
 *             ---->  
 *  y:  k  l             y:  j  l
 *
 * Suppose, for example, that x=2, or in binary, x=0b10. Furthermore,
 * suppose that y=3=0b11. Then after a call to trBits, we'll have
 * x=3=0b11 and y=1=0b01.
 *
 * Any higher-order bits are ignored. For example, the preceding
 * example still generates x=3 and y=1, even if on input y=7=0b111.
 */
static inline void trBits(unsigned char *x, unsigned char *y);

static inline void trBits(unsigned char *x, unsigned char *y) {
    unsigned char tmp;

    tmp = (*x & 2) | ((*y & 2) >> 1);
    *y = ((*x & 1) << 1) | (*y & 1);
    *x = tmp;
}

static inline unsigned char gamete(unsigned char x, unsigned char y);
static inline unsigned char gamete(unsigned char x, unsigned char y) {
    assert((x & ~1u) == '\0');  /* x == 0 or 1 */
    assert((y & ~1u) == '\0');  /* y == 0 or 1 */
    return ((x << 1) | y);
}

#endif
