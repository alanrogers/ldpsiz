/**
 * @file snp.c
 * @author Alan R. Rogers
 * @brief A single nucleotide polymorphism (SNP).
 * 
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "snp.h"
#include "misc.h"
#include "sums.h"
#include "em.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/**
 * SNPstore is a place to store objects of type SNP when they are not
 * in use. This allows them to be reused without repeated calls to
 * malloc.    
 */
struct SNPstore {
    unsigned    nGtype;/**< number of genotypes at each SNP */
    int         bootreps;
                       /**< number of bootstrap replicates */
    SNP        *head;  /**< pointer to first SNP */
};

double      SNP_getDsq_haploid(double *pqpq, SNP * x, SNP * y);

/**
 * Set values to initial states. This function does not initialize
 * snp->prev, because doing so risks a memory leak.
 */
void SNP_clear(SNP * snp) {
    if(snp == NULL)
        return;
    snp->mappos = -1.0;
    snp->ndx = -1;
    snp->p = -1.0;
    snp->pq = -1.0;
    assert(snp->prev == NULL);
    return;
}

/** free linked list of SNPs */
void SNP_free(SNP * snp) {
    if(snp == NULL)
        return;
    SNP_free(snp->prev);
#ifndef NDEBUG
    if(snp->multiplicity)
        myassert(snp->bootreps > 0);
#endif
    free(snp);
    return;
}

/**
 * Allocate a new SNP, including the arrays used to store
 * genotype data and bootstrap multiplicity values. In order to
 * improve efficiency, these allocations are all done with a single
 * call to malloc. That way, data for each SNP is in contiguous
 * memory, and we minimize cache misses.
 *
 * @param [in] nGtype Number of genotypes at each SNP
 * @param [in] bootreps Number of bootstrap replicates.
 * @returns pointer to newly-allocated SNP.
 */
SNP        *SNP_new(unsigned nGtype, int bootreps) {
    SNP        *snp;
    static const SNP dummy;     /* used only in sizeof arguments */

    size_t      size = sizeof(SNP);

    size += nGtype * sizeof(dummy.gtype[0]);
    if(bootreps > 0)
        size += bootreps * sizeof(dummy.multiplicity[0]);

    snp = malloc(size);
    checkmem(snp, __FILE__, __LINE__);

    snp->gtype = (void *) &(snp[1]);

    if(bootreps > 0)
        snp->multiplicity = (void *) &(snp->gtype[nGtype]);
    else
        snp->multiplicity = NULL;

    snp->nGtype = nGtype;
    snp->bootreps = bootreps;
    snp->prev = NULL;

    SNP_clear(snp);

#ifndef NDEBUG
    size_t      tst;

    if(bootreps)
        tst = (size_t) & (snp->multiplicity[bootreps]);
    else
        tst = (size_t) & (snp->gtype[nGtype]);
    tst -= (size_t) snp;
    /* "tst" should equal "size" */
    if(size != tst) {
        fprintf(stderr, "\nERR@%s:%d in proc SNP_new:"
                " size=%lu != tst=%lu\n\n",
                __FILE__, __LINE__,
                (long unsigned) size, (long unsigned) tst);
        exit(EXIT_FAILURE);
    }
#endif

    return (snp);
}

/**
 * Initialize a SNP using the data given by the other parameters.
 *
 * @param [in] ndx SNPs are indexed (numbered) starting with 0. ndx is the
 * index of the current SNP.
 *
 * @param [in] mappos position of SNP in units of the recombinational
 * map.
 *
 * @param [in] gtype points to a character string representing
 * genotypes
 *
 * @param [in,out] boot points to Boot structure, which is used to
 * store information about bootstrap replicates.
 *
 * @param [in] ploidy (either 0 or 1)
 *
 * @returns 1 if the SNP is polymorphic, or 0 if monomorphic.
 */
int SNP_set(SNP * snp, long ndx, double mappos,
            const unsigned char *gtype, const Boot * boot, unsigned ploidy) {
    long        i;

    snp->mappos = mappos;
    snp->ndx = ndx;
    if(ploidy == 2)
        snp->sum = sumDiploid(gtype, snp->nGtype);
    else {
        assert(ploidy == 1);
        snp->sum = sum_char(gtype, snp->nGtype);
    }
    snp->p = snp->sum / ((double) (ploidy * snp->nGtype));

    /*
     * Perhaps I should multiply pq by n/(n-1) to correct for bias.
     */
    snp->pq = snp->p * (1.0 - snp->p);

    memcpy(snp->gtype, gtype, snp->nGtype * sizeof(snp->gtype[0]));

    if(snp->bootreps > 0) {
        myassert(snp->multiplicity);
        myassert(boot);
        myassert(Boot_nReps(boot) == snp->bootreps);
        for(i = 0; i < snp->bootreps; ++i)
            snp->multiplicity[i] = Boot_multiplicity(boot, ndx, i);
    } else {
        myassert(snp->multiplicity == NULL);
        myassert(boot == NULL);
    }
    if(snp->pq > 0.0)
        return 1;
    return 0;
}

/** Print SNP */
void SNP_show(SNP * snp, FILE * fp) {
    int         i;

    fprintf(fp, "SNP at %lx:\n", (unsigned long) snp);
    fprintf(fp, "    mappos: %lg\n", snp->mappos);
    fprintf(fp, "    nGtype: %d\n", snp->nGtype);
    fprintf(fp, "    ndx   : %ld\n", snp->ndx);
    fprintf(fp, "    p     : %g\n", snp->p);
    fprintf(fp, "    pq    : %g\n", snp->pq);
    fprintf(fp, "    prev  : %lx\n", (unsigned long) snp->prev);
    fprintf(fp, "    ");
    for(i = 0; i < snp->nGtype; ++i)
        fprintf(fp, " %1u", snp->gtype[i]);
    putc('\n', fp);
}

/** Count SNPs in linked list */
int SNP_count(SNP * snp) {
    if(snp == NULL)
        return 0;
    return 1 + SNP_count(snp->prev);
}

/* Return map position of SNP */
double SNP_mappos(const SNP * snp) {
    assert(snp);
    return snp->mappos;
}

/**
 * Connect list2 to the tail of list1 and return a pointer to list1.
 *
 * Usage: list1 = SNP_connect(list1, list2)
 *
 * @param[in] list1,list2 Linked lists of SNPs, which will be
 * concatenated.
 *
 * @returns A pointer to the head of the combined list. If list1 is
 * nonNULL, this returned pointer will equal "list1". Otherwise, it will
 * equal "list2". Either way, the returned value is the head of a list
 * that contains all the entries in list 1 followed by all the entries
 * in list 2.
 */
SNP        *SNP_connect(SNP * list1, SNP * list2) {
    SNP        *snp;

    if(list1 == NULL)
        return list2;

    for(snp = list1; snp->prev != NULL; snp = snp->prev)
        ;

    myassert(snp != NULL);
    myassert(snp->prev == NULL);

    snp->prev = list2;
    return list1;
}

/// Count copies of allele with value 1
unsigned SNP_countDerived(SNP * snp, unsigned ploidy) {
    if(ploidy == 2)
        return sumDiploid(snp->gtype, snp->nGtype);
    return sum_char(snp->gtype, snp->nGtype);
}

/// Count copies of minor allele.
unsigned SNP_countMinor(SNP * snp, unsigned ploidy) {
    unsigned count = SNP_countDerived(snp, ploidy);
    unsigned complement = ploidy*snp->nGtype - count;
    return (count<complement ? count : complement);
}

/**
 * Calculate D squared from pair of haploid SNPs.
 *
 * Use two SNPs to calculate D*D (where D is the standard LD
 * coefficient). Also calculates the product of the pq values at the
 * two SNPs. Function returns the value of D*D and sets the value of
 * *pqpq. This version is for haploid genotypes, with each genotype
 * equal either to 0 or 1.
 */
double SNP_getDsq_haploid(double *pqpq, SNP * x, SNP * y) {
    double      D = 0.0;

    /* calculate sum(x[i] * y[i]) */
    D = sum_and_char(x->gtype, y->gtype, x->nGtype);

    D -= x->sum * y->p;
    D /= (x->nGtype - 1);

    *pqpq = x->pq * y->pq;

    return D * D;
}

/**
 * Calculate D squared.
 *
 * Use two SNPs to calculate D*D (where D is the standard LD
 * coefficient) and the product of the pq values at the two
 * SNPs. Function also sets the value of *pqpq.  This version uses the
 * EM algorithm to get a maximum-likelihood estimate of D. This also
 * provides a ML estimate of Dsq=D*D because of the invariant property
 * of ML estimates.
 */
double SNP_getDsq(double *pqpq, SNP * x, SNP * y, unsigned ploidy) {

    myassert(x->nGtype == y->nGtype);
    int         status;
    register unsigned j = 0;

    *pqpq = x->pq * y->pq;
    myassert(*pqpq > 0.0);

    DsqData     dd = {
        .tol = sqrt(DBL_EPSILON),
        .alpha = *pqpq,
        .beta = x->p + y->p - 2.0 * x->p * y->p,
        .px = x->p,
        .py = y->p,
        .nGtype = x->nGtype,
        .nUnphased = 0,
        .nGam = {0}
    };

    /*
     * Tabulate gametes into nGam array. Haploid version has unrolled
     * loops for speed. Diploid version also counts the number of
     * unphased double heterozygotes.
     */
    if(ploidy == 2) {
        for(j = 0; j < x->nGtype; ++j) {
            if(UNPHASED_PAIR(x->gtype[j], y->gtype[j]))
                ++dd.nUnphased;
            else {
                /* Count phased gametes of each type */
                unsigned char g1 = MAKE_PHASED(x->gtype[j]);
                unsigned char g2 = MAKE_PHASED(y->gtype[j]);

                trBits(&g1, &g2);
                ++dd.nGam[g1];
                ++dd.nGam[g2];
            }
        }
    } else {
        /* haploid */
        /* unrolled loop */
        register unsigned k = x->nGtype % 5;

        for(j = 0; j < k; ++j)
            ++dd.nGam[gamete(x->gtype[j], y->gtype[j])];

        for(j = k; j < x->nGtype; j += 5) {
            ++dd.nGam[gamete(x->gtype[j], y->gtype[j])];
            ++dd.nGam[gamete(x->gtype[j + 1], y->gtype[j + 1])];
            ++dd.nGam[gamete(x->gtype[j + 2], y->gtype[j + 2])];
            ++dd.nGam[gamete(x->gtype[j + 3], y->gtype[j + 3])];
            ++dd.nGam[gamete(x->gtype[j + 4], y->gtype[j + 4])];
        }
    }

    dd.loD = loD(dd.px, dd.py, dd.nGam);
    dd.hiD = hiD(dd.px, dd.py, dd.nGam);

    double      D;

    status = minimize1D(&D, &dd);
    if(status) {
        DsqData_print(&dd, __FILE__, __LINE__, stdout);
        eprintf("ERR@%s:%d: minimize1D returned %d\n",
                __FILE__, __LINE__, status);
    }

    return D * D;
}

/** index of this SNP (begins with 0) */
long SNP_ndx(const SNP * snp) {
    myassert(snp);
    return snp->ndx;
}

/** Number of genotypes sampled for this SNP */
int SNP_nGtype(const SNP * snp) {
    myassert(snp);
    return snp->nGtype;
}

/**
 * Allocate a new object of type SNPstore.
 *
 * @param nGtype number of genotypes
 * @param bootreps number of bootstrap replicates
 */
SNPstore   *SNPstore_new(unsigned nGtype, int bootreps) {
    SNPstore   *store = (SNPstore *) malloc(sizeof(SNPstore));

    checkmem(store, __FILE__, __LINE__);
    store->nGtype = nGtype;
    store->bootreps = bootreps;
    store->head = NULL;
    return (store);
}

/** Free an object of type SNPstore */
void SNPstore_free(SNPstore * store) {
    SNP_free(store->head);
    free(store);
    return;
}

/** Obtain a SNP from the store, allocating if necessary. */
SNP        *SNPstore_checkout(SNPstore * store) {
    SNP        *snp;

    if(store->head == NULL)
        return SNP_new(store->nGtype, store->bootreps);
    snp = store->head;
    store->head = snp->prev;
    snp->prev = NULL;
    SNP_clear(snp);
    return (snp);
}

/** Check a SNP back into the store. */
void SNPstore_checkin(SNPstore * store, SNP * snp) {
    store->head = SNP_connect(snp, store->head);
    return;
}

#ifndef NDEBUG
#include <time.h>
void SNP_test(int verbose) {

    unsigned    nGtype = 6;
    int         rval;
    unsigned char gtype1[] = { 0, 0, 0, 1, 1, 1 };
    unsigned char gtype2[] = { 0, 0, 1, 1, 1, 1 };
    int         bootreps = 0;
    Boot       *boot = NULL;
    const int   folded = true;

    SNP        *snp1 = SNP_new(nGtype, bootreps);

    assert(SNP_nGtype(snp1) == nGtype);
    assert(snp1->bootreps == bootreps);
    assert(1 == SNP_count(snp1));
    assert(-1.0 == SNP_mappos(snp1));
    assert(-1 == snp1->ndx);
    assert(-1.0 == snp1->p);
    assert(-1.0 == snp1->pq);
    if(bootreps > 0)
        assert(snp1->multiplicity);
    else
        assert(NULL == snp1->multiplicity);
    assert(NULL == snp1->prev);

    if(verbose)
        unitTstResult("SNP_new", "OK");

    long        ndx = 0;
    double      mappos = 0.005;
    unsigned    ploidy = 1;
    unsigned    twoNsmp = nGtype * ploidy;

    rval = SNP_set(snp1, ndx, mappos, gtype1, boot, ploidy);
    assert(rval == 1);
    assert(SNP_ndx(snp1) == ndx);
    assert(SNP_mappos(snp1) == mappos);
    assert(0.5 == snp1->p);
    assert(0.25 == snp1->pq);
    assert(3 == snp1->sum);

    if(verbose)
        SNP_show(snp1, stdout);

    SNP        *snp2 = SNP_new(nGtype, bootreps);

    ndx = 1;
    mappos = 0.01;
    rval = SNP_set(snp2, ndx, mappos, gtype2, boot, ploidy);
    assert(rval == 1);
    assert(SNP_ndx(snp2) == ndx);
    assert(SNP_mappos(snp2) == mappos);

    if(verbose)
        SNP_show(snp2, stdout);

    double      tol = sqrt(DBL_EPSILON);
    double      Dsq, pqpq;
    double      EDsq = 0.027777777777777783;
    double      Epqpq = 0.05555555555555556;

    Dsq = SNP_getDsq(&pqpq, snp1, snp2, ploidy);
    if(verbose)
        printf("haploid: Dsq=%lg (expected %lg); pqpq=%lg (expected %lg)\n",
               Dsq, EDsq, pqpq, Epqpq);
    assert(fabs(Dsq - EDsq) < tol);
    assert(fabs(pqpq - Epqpq) < tol);

    if(verbose)
        unitTstResult("SNP_getDsq (haploid)", "OK");

    /* test diploid case */
    ploidy = 2;
    ndx = 0;
    mappos = 0.005;
    unsigned char gtype3[] = { 0, 0, 1, 3, 3, 3 };
    unsigned char gtype4[] = { 0, 0, 1, 3, 3, 3 };
    rval = SNP_set(snp1, ndx, mappos, gtype3, boot, ploidy);
    ndx = 1;
    mappos = 0.01;
    rval = SNP_set(snp2, ndx, mappos, gtype4, boot, ploidy);

    double      EDsq1 = 0.05907600308641974;

    Epqpq = 0.05907600308641975;
    Dsq = SNP_getDsq(&pqpq, snp1, snp2, ploidy);
    if(verbose)
        printf("diploid1: Dsq=%lg (expected %lg); pqpq=%lg (expected %lg)\n",
               Dsq, EDsq1, pqpq, Epqpq);
    assert(fabs(Dsq - EDsq1) < tol);
    assert(fabs(pqpq - Epqpq) < tol);

    unsigned char gtype5[] = { 0, 0, 2, 3, 3, 3 };
    rval = SNP_set(snp1, ndx, mappos, gtype5, boot, ploidy);

    double      EDsq2 = 0.025511188271604916;

    Dsq = SNP_getDsq(&pqpq, snp1, snp2, ploidy);
    if(verbose)
        printf("diploid2: Dsq=%lg (expected %lg); pqpq=%lg (expected %lg)\n",
               Dsq, EDsq2, pqpq, Epqpq);
    assert(fabs(Dsq - EDsq2) < tol);
    assert(fabs(pqpq - Epqpq) < tol);

    unsigned char gtype6[] = { 0, 0, 4, 3, 3, 3 };
    rval = SNP_set(snp1, ndx, mappos, gtype6, boot, ploidy);
    Dsq = SNP_getDsq(&pqpq, snp1, snp2, ploidy);
    if(verbose)
        printf("diploid3: Dsq=%lg (expected %lg-%lg);"
               " pqpq=%lg (expected %lg)\n", Dsq, EDsq2, EDsq1, pqpq, Epqpq);
    assert(Dsq >= EDsq2);
    assert(Dsq <= EDsq1);
    assert(fabs(pqpq - Epqpq) < tol);

    if(verbose)
        unitTstResult("SNP_getDsq (diploid)", "OK");

    bootreps = 10;
    long        nSNPs = 100;
    long        blockLength = 1;
    double      windowcm = 0.2;
    int         nbins = 10;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

    gsl_rng_set(rng, (unsigned) time(NULL));
    boot = Boot_new(nSNPs, bootreps, twoNsmp, folded, blockLength,
                    windowcm, nbins, rng);
    SNP        *snp3 = SNP_new(nGtype, bootreps);

    rval = SNP_set(snp3, ndx, mappos, gtype1, boot, ploidy);
    long        nBlocks = round(nSNPs / ((double) blockLength));
    unsigned    i;

    for(i = 0; i < bootreps; ++i) {
        long        m = SNP_multiplicity(snp3, i);

        if(verbose)
            printf("multiplicity[%2u]=%ld\n", i, m);
        assert(m <= nBlocks);
        assert(0 <= m);
    }

    SNP_free(snp1);
    SNP_free(snp2);
    SNP_free(snp3);
    Boot_free(boot);
    snp1 = snp2 = snp3 = NULL;
    boot = NULL;

    unitTstResult("SNP", "OK");

    ploidy = 1;
    bootreps = 0;
    snp1 = SNP_new(nGtype, bootreps);
    snp2 = SNP_new(nGtype, bootreps);

    SNPstore   *store = SNPstore_new(nGtype, bootreps);

    assert(store);
    assert(nGtype == store->nGtype);
    assert(bootreps == store->bootreps);
    assert(NULL == store->head);

    if(verbose)
        unitTstResult("SNPstore_new", "OK");

    SNPstore_checkin(store, snp1);
    snp1 = NULL;
    assert(1 == SNP_count(store->head));
    SNPstore_checkin(store, snp2);
    snp2 = NULL;
    assert(2 == SNP_count(store->head));

    if(verbose) {
        unitTstResult("SNP_count", "OK");
        unitTstResult("SNPstore_checkin", "OK");
    }

    snp1 = SNPstore_checkout(store);
    assert(nGtype == SNP_nGtype(snp1));
    assert(bootreps == snp1->bootreps);
    assert(1 == SNP_count(store->head));

    SNPstore_free(store);
    SNP_free(snp1);
    snp1 = NULL;
    store = NULL;

    if(verbose)
        unitTstResult("SNPstore_free", "OK");

    unitTstResult("SNPstore", "OK");

    gsl_rng_free(rng);
}
#endif
