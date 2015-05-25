/**
 * @file snp.h
 * @author Alan R. Rogers
 * @brief Header for snp.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SNP_H
#define LDPSIZ_SNP_H

#include <stdio.h>
#include "typedefs.h"
#include "boot.h"
#include "em.h"
#include "misc.h"

/**
 * The SNP class represents a single nucleotide polymorphism.
 *
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
struct SNP {
    double      mappos;         /* position of SNP in centimorgans */
    unsigned char *gtype;       /* binary-encoded genotypes in sample */
    unsigned    nGtype;         /* number of entries in array gtype */
    long        ndx;            /* index of current SNP (begins with 0) */
    double      sum, p, pq;
    int         bootreps;
    unsigned   *multiplicity;
    struct SNP *prev;           /* previous SNP in linked list */
};

void        SNP_clear(SNP * snp);
SNP        *SNP_connect(SNP * list1, SNP * list2);
int         SNP_count(SNP * snp);
void        SNP_free(SNP * snp);
unsigned    SNP_countDerived(SNP * snp, unsigned ploidy);
unsigned    SNP_countMinor(SNP * snp, unsigned ploidy);
double      SNP_getDsq(double *pqpq, SNP * x, SNP * y, unsigned ploidy);
double      SNP_mappos(const SNP * snp);
static      inline unsigned SNP_multiplicity(const SNP * snp, int rep);
long        SNP_ndx(const SNP * snp);
SNP        *SNP_new(unsigned nGtype, int bootreps);
int         SNP_nGtype(const SNP * snp);
int         SNP_set(SNP * snp, long ndx, double mappos,
                    const unsigned char *gtype, const Boot * boot,
                    unsigned ploidy);
void        SNP_show(SNP * snp, FILE * fp);

void        SNPstore_checkin(SNPstore * store, SNP * snp);
SNP        *SNPstore_checkout(SNPstore * store);
void        SNPstore_free(SNPstore * sp);
SNPstore   *SNPstore_new(unsigned nGtype, int bootreps);

/* inline functions */
static      inline unsigned SNP_multiplicity(const SNP * snp, int rep) {
    myassert(snp);
    myassert(rep < snp->bootreps);
    return snp->multiplicity[rep];
}

#  ifndef NDEBUG
void        SNP_test(int verbose);
#  endif
#endif
