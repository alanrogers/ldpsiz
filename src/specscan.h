/**
 * @file specscan.h
 * @author Alan R. Rogers
 * @brief Header for specscan.c
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SPECSCAN_H
#define LDPSIZ_SPECSCAN_H

#include <stdio.h>
#include "typedefs.h"
#include "boot.h"
#include "em.h"
#include "misc.h"

/**
 * Estimate the site frequency spectrum. Definition is provided openly
 * here (rather than hidden in specscan.c) so that functions can be
 * inlined.  
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
    int        *multiplicity;
    struct SNP *prev;           /* previous SNP in linked list */
};

void        SNP_clear(SNP * snp);
SNP        *SNP_connect(SNP * list1, SNP * list2);
int         SNP_count(SNP * snp);
void        SNP_free(SNP * snp);
double      SNP_getSpec(unsigned dim, double *spec, SNP * x, unsigned ploidy);
double      SNP_mappos(const SNP * snp);
static      inline int SNP_multiplicity(const SNP * snp, int rep);
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

int         Specscan_advance(Specscan * specscan, Tabulation * tab, Boot * boot,
                           long lineno);
SNP        *Specscan_currSNP(Specscan * specscan);
void        Specscan_free(Specscan * specscan);
int         Specscan_nextSNP(Specscan * specscan, Boot * boot);
Specscan     *Specscan_new(double width_cm, FILE * ifp,
                       long sampling_interval, unsigned ploidy);
unsigned    Specscan_nGtype(const Specscan * specscan);
long        Specscan_nSNPsRead(const Specscan * specscan);

#ifndef NDEBUG
void        Specscan_test(int verbose);
#endif

/* inline functions */
static      inline int SNP_multiplicity(const SNP * snp, int rep) {
    myassert(snp);
    myassert(rep < snp->bootreps);
    return snp->multiplicity[rep];
}

#endif
