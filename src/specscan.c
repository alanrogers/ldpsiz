/**
 * @file specscan.c
 * @author Alan R. Rogers
 * @brief A window that slides across the chromosome.
 * 
 * LD is calculated from pairs of sites within a window, which
 * slides across the chromosome. This file implements Window, which
 * represents that sliding window.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include "specscan.h"
#include "readgtp.h"
#include "misc.h"
#include "tabulation.h"
#include "boot.h"
#include "sums.h"

/**
 * @brief Allocate a new Window.
 *
 * @param [in] width_cm width of Window in centimorgans
 *
 * @param [in] sampling_interval If sampling_interval is 2, window
 * will skip every other SNP. This makes things fast for debugging.
 *
 * @param [in] ploidy 1 for haploid, 2 for diploid
 */
Window     *Window_new(double width_cm,
                       FILE * ifp,
                       long int sampling_interval, unsigned ploidy) {
    Window     *window = (Window *) malloc(sizeof(Window));

    checkmem(window, __FILE__, __LINE__);

    window->ifp = ifp;
    window->width_cm = width_cm;
    window->curr = NULL;
    window->nSNPs = 0;
    window->nGtype = 0;
    window->ploidy = ploidy;
    window->store = NULL;
    window->sampling_interval = sampling_interval;
    return (window);
}

/** Free a Window */
void Window_free(Window * window) {
    if(window->curr)
        SNP_free(window->curr);
    if(window->store)
        SNPstore_free(window->store);
    free(window);
    return;
}

/** Return pointer to current focal SNP */
SNP        *Window_currSNP(Window * window) {
    assert(window);
    return window->curr;
}

/**
 * Get next SNP from input file and link it into the linked list
 * pointed to by window->curr. Return 0 on success, EOF if end of file
 * is reached, and 1 on any other error.
 */
int Window_nextSNP(Window * window, Boot * boot) {
    double      mappos;
    int         nGtype, polymorphic = 0;
    SNP        *snp = NULL;

    /* get next SNP from file, skipping duplicate mappos values */
    do {
        nGtype = Gtp_readSNP(window->ifp, NULL, 0,  /* don't get snpId */
                             &mappos, NULL, 0,  /* don't get list of alleles */
                             window->gtype, sizeof(window->gtype),
                             window->ploidy == 2);
    } while(nGtype != EOF && window->curr
            && Dbl_near(mappos, window->curr->mappos));
    if(nGtype == EOF)
        return EOF;

    if(window->nGtype == 0) {
        /* Initialization code only runs once per thread. */
        window->nGtype = nGtype;
        window->store = SNPstore_new(window->nGtype,
                                     (boot ? Boot_nReps(boot) : 0));
    }
#ifndef NDEBUG
    if(window->nGtype != nGtype) {
        fprintf(stderr, "ERR@%s:%d: window->nGtype(=%d) != nGtype(=%u))\n",
                __FILE__, __LINE__, window->nGtype, nGtype);
        dostacktrace(__FILE__, __LINE__, stderr);
        exit(1);
    }
#endif

    if(snp == NULL)
        snp = SNPstore_checkout(window->store);

    polymorphic = SNP_set(snp, window->nSNPs, mappos,
                          window->gtype, boot, window->ploidy);
    if(!polymorphic) {
		fflush(stdout);
        fprintf(stderr, "%s:%d: Gtp_readSNP returned a monomorphic SNP",
                __FILE__, __LINE__);
		SNP_show(snp, stderr);
		exit(1);
	}

    /* link new SNP into list */
    snp->prev = window->curr;
    window->curr = snp;
    ++window->nSNPs;

    return 0;
}

/** Add a new SNP to the window. */
int Window_advance(Window * window, Tabulation * tab, Boot * boot, long count) {
    double      sep_cm;
    int         rval;
    SNP        *snp;

    rval = Window_nextSNP(window, boot);
    if(rval == EOF)
        return EOF;

    /*
     * If sampling_interval > 1, the following code will seldom
     * execute. This increases speed by roughly a factor of
     * sampling_interval.
     */
    if(count % window->sampling_interval == 0) {
        /* Each pass through loop compares current SNP with a previous
         * SNP, provided that previous SNP is within the window. If it is
         * outside the window, the list is truncated at that point, and
         * SNPs outside the window are returned to the store.
         */
        for(snp = window->curr; snp->prev != NULL; snp = snp->prev) {

            /* ignore pairs with zero separation */
            if(Dbl_near(window->curr->mappos, snp->prev->mappos))
                continue;

            sep_cm = window->curr->mappos - snp->prev->mappos;
            myassert(sep_cm > 0.0);

            /* If we've reached the window, then truncate the list */
            if(sep_cm >= window->width_cm) {
                SNPstore_checkin(window->store, snp->prev);
                snp->prev = NULL;
                break;
            }

            double      Dsq, pqpq;

            /* Dsq and pqpq are values to be tabulated */
            Dsq = SNP_getDsq(&pqpq, window->curr, snp->prev, window->ploidy);

            /* record these values with weight 1 */
            Tabulation_record(tab, Dsq, pqpq, sep_cm, 1);

            if(boot)
                Boot_addLD(boot, Dsq, pqpq, sep_cm, window->curr, snp->prev);
        }
    }
    return 0;
}

/** Return current value of window->nGtype */
unsigned Window_nGtype(const Window * window) {
    return window->nGtype;
}

/** Return the number of SNPs that have been read. */
long Window_nSNPsRead(const Window * window) {
    myassert(window);
    return window->nSNPs;
}

#ifndef NDEBUG
void Window_test(int verbose) {
    const char *tstInput = "# source              = test input\n\
# Haploid sample size   = 10\n\
# Ploidy                = %d\n\
#   snp_id     nucpos    mappos alleles genotypes\n\
         0        262    0.0262       01 0001\n\
         1        362    0.0362       01 1100\n\
         2        536    0.0536       01 1000\n\
         3        799    0.0799       01 0010\n\
         4        861    0.0861       01 0010\n\
         5       1337    0.1337       01 0110\n\
         6       1564    0.1564       01 1110\n\
         7       1905    0.1905       01 0010\n\
         8       1968    0.1968       01 1001\n\
         9       2419    0.2419       01 0010\n";

    unsigned    nGtype = 6;
    int         rval;
    unsigned char gtype1[] = { 0, 0, 0, 1, 1, 1 };
    unsigned char gtype2[] = { 0, 0, 1, 1, 1, 1 };
    int         bootreps = 0;

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

    Boot       *boot = NULL;

    long        ndx = 0;
    double      mappos = 0.005;
    unsigned    ploidy = 1;

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
    boot = Boot_new(nSNPs, bootreps, blockLength, windowcm, nbins, rng);
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

    const char *fname = "window-tmp.gtp";
    FILE       *fp = fopen(fname, "w");

    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(fname, "r");

    windowcm = 0.1;
    long        sampling_interval = 1;

    ploidy = 1;
    Window     *window = Window_new(windowcm, fp, sampling_interval, ploidy);

    assert(0 == Window_nGtype(window));
    assert(0 == Window_nSNPsRead(window));
    assert(windowcm == window->width_cm);
    assert(ploidy == window->ploidy);
    assert(NULL == Window_currSNP(window));

    /* advance past header material */
    Assignment *asnmt = Gtp_readHdr(fp);

    Assignment_free(asnmt);
    asnmt = NULL;

    long        lineno = 0;

    Window_nextSNP(window, boot);

    assert(4 == Window_nGtype(window));
    assert(1 == Window_nSNPsRead(window));
    assert(NULL != window->store);
    assert(NULL != Window_currSNP(window));

    Tabulation *tab = Tabulation_new(windowcm, nbins);

    ++lineno;
    Window_advance(window, tab, boot, lineno);

    assert(4 == Window_nGtype(window));
    assert(2 == Window_nSNPsRead(window));
    assert(NULL != window->store);
    assert(NULL != Window_currSNP(window));

    snp1 = Window_currSNP(window);

    assert(Dbl_near(0.0362, snp1->mappos));
    assert(1 == snp1->ndx);
    assert(4 == snp1->nGtype);
    assert(1 == snp1->gtype[0]);
    assert(1 == snp1->gtype[1]);
    assert(0 == snp1->gtype[2]);
    assert(0 == snp1->gtype[3]);
    assert(Dbl_near(0.0262, snp1->prev->mappos));

    Window_free(window);
    Tabulation_free(tab);
    tab = NULL;
    window = NULL;

    unitTstResult("Window", "OK");

    gsl_rng_free(rng);
    fclose(fp);
}
#endif
