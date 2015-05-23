/**
 * @file window.c
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

#include "window.h"
#include "readgtp.h"
#include "misc.h"
#include "tabulation.h"
#include "spectab.h"
#include "boot.h"
#include "sums.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

/**
 * Window represents a window that slides across the chromosome.
 */
struct Window {
    FILE       *ifp;            /**< input file pointer */
    double      width_cm;       /**< the width of the window in cM */
    SNP        *curr;           /**< points to current focal SNP */
    long        nSNPs;          /**< number of SNPs that have been read */
    unsigned    ploidy;         /**< 1 for haploid, 2 for diploid */
    unsigned char gtype[500];   /**< binary-encoded genotypes in sample */
    unsigned    nGtype;         /**< number of genotypes in data */
    SNPstore   *store;          /**< holds SNPs when not in use */
    long int    sampling_interval;
                                /**< for debugging */
};

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
        // Initialization code only runs once per thread.
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

/// Add a new SNP to the window.
int Window_advance(Window * window, Tabulation * tab, Spectab *spectab,
                   Boot * boot, long lineno) {
    double      sep_cm;
    int         rval;
    SNP        *snp;

    rval = Window_nextSNP(window, boot);
    if(rval == EOF)
        return EOF;

    // If sampling_interval > 1, the following code will seldom
    // execute. This increases speed by roughly a factor of
    // sampling_interval.
    if(lineno % window->sampling_interval == 0) {
        // Add current SNP to site frequency spectrum.
        // weight=1 because this isn't a bootstrap replicate.
        unsigned    alleleCount = SNP_countMinor(window->curr,
                                                 window->ploidy);
        Spectab_record(spectab, alleleCount, 1);

        // Each pass through loop compares current SNP with a previous
        // SNP, provided that previous SNP is within the window. If it is
        // outside the window, the list is truncated at that point, and
        // SNPs outside the window are returned to the store.
        for(snp = window->curr; snp->prev != NULL; snp = snp->prev) {

            // ignore pairs with zero separation
            if(Dbl_near(window->curr->mappos, snp->prev->mappos))
                continue;

            sep_cm = window->curr->mappos - snp->prev->mappos;
            myassert(sep_cm > 0.0);

            // If we've reached the window, then truncate the list
            if(sep_cm >= window->width_cm) {
                SNPstore_checkin(window->store, snp->prev);
                snp->prev = NULL;
                break;
            }

            double      Dsq, pqpq;

            // Dsq and pqpq are values to be tabulated
            Dsq = SNP_getDsq(&pqpq, window->curr, snp->prev, window->ploidy);

            // record these values with weight 1
            Tabulation_record(tab, Dsq, pqpq, sep_cm, 1);

            if(boot)
                Boot_addLD(boot, Dsq, pqpq, sep_cm, window->curr, snp->prev);
        }
    }
    return 0;
}

/// Return window->nGtype
unsigned Window_nGtype(const Window * window) {
    return window->nGtype;
}

/// Return window->ploidy
unsigned    Window_ploidy(const Window * window) {
    return window->ploidy;
}


/// Return the number of SNPs that have been read.
long Window_nSNPsRead(const Window * window) {
    myassert(window);
    return window->nSNPs;
}

#ifndef NDEBUG
#include "assign.h"
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

    const char *fname = "window-tmp.gtp";
    FILE       *fp = fopen(fname, "w");

    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(fname, "r");

    long        blockLength = 1;
    long        nSNPs = 100;
    double      windowcm = 0.1;
    int         nbins = 10;
    long        sampling_interval = 1;
    unsigned    ploidy;
    int         folded = true;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned) time(NULL));
    Boot       *boot = Boot_new(nSNPs, 0 /* bootreps=0 */, blockLength,
                                windowcm, nbins, rng);
    assert(boot == NULL);

    int         bootreps = 10;
    boot = Boot_new(nSNPs, bootreps, blockLength, windowcm, nbins, rng);
    assert(boot != NULL);

    ploidy = 1;
    Window     *window = Window_new(windowcm, fp, sampling_interval, ploidy);

    assert(0 == Window_nGtype(window));
    assert(ploidy = Window_ploidy(window));
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
    unsigned nHapSamp = Window_nGtype(window) * Window_ploidy(window);
    Spectab *spectab = Spectab_new(nHapSamp, folded);

    ++lineno;
    Window_advance(window, tab, spectab, boot, lineno);

    Spectab_print(spectab, stdout);

    assert(4 == Window_nGtype(window));
    assert(2 == Window_nSNPsRead(window));
    assert(NULL != window->store);
    assert(NULL != Window_currSNP(window));

    SNP        *snp1 = Window_currSNP(window);

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
    Spectab_free(spectab);
    tab = NULL;
    window = NULL;
    spectab = NULL;

    unitTstResult("Window", "OK");

    gsl_rng_free(rng);
    fclose(fp);
}
#endif
