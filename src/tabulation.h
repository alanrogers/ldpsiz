/**
 * @file tabulation.h
 * @author Alan R. Rogers
 * @brief Header for tabulation.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_TABULATION_H
#define LDPSIZ_TABULATION_H

#include "typedefs.h"
#include "misc.h"
#include "array.h"

/*
 * tabulation.h. This file defines the Tabulation object, which keeps
 * track of comparisons between pairs of SNPs. These comparisons take
 * place within a window of size windocm centimorgans, which slides
 * across the genome. Comparisons are tabulated into nbins bins, based
 * on the distance that separates them along the chromosome. The
 * Tabulation object represents a table, with a row for each bin, and
 * a column for each of: nobs, sep_cm, numerator, denominator. Each
 * time a pair of SNPs is compared, we add something to each of these
 * columns in the row determined by the distance that separates the
 * SNPs along the chromosome.
 *
 * The definition of struct Tabulation is exposed here (rather than
 * hidden in tabulation.c) for the benefit of the inline functions
 * defined at the bottom of this file.
 */
struct Tabulation {
    int         nbins;          /* how many bins to use in tabulating results */
    double      windowcm;       /* width of sliding window in centimorgans */
    double      bins_per_cm;    /* nbins / windowcm */
    long unsigned *nobs;        /* number of observations w/i each bin */
    double     *sep_cm;         /* summed separation of SNPs, by bin */
    double     *numerator;      /* summed D^2 values */
    double     *denominator;    /* summed values of p_A q_A p_B q_B */
    double     *rsq;            /* summed values of r^2 */

    /*
     * Hold previous values of nobs[0] and nobs[nbins-1]. If new value
     * is less than previous value, then an overflow has occurred.
     */
    long unsigned oldNobs[2];
};

Tabulation *Tabulation_new(double windowcm, int nbins);
Tabulation *Tabulation_dup(Tabulation * old);
void        Tabulation_free(Tabulation * tab);
void        Tabulation_plus_equals(Tabulation * x, const Tabulation * y);
void        Tabulation_print(Tabulation * tab, FILE * ofp);
#if 0
int         Tabulation_report(Tabulation * tab, double *sep_cm, long unsigned *nobs,
                              double *sigdsq, double *rsq);
#else
int         Tabulation_report(Tabulation *tab,
                              DblArray *sep_cm, ULIntArray *nobs,
                              DblArray *sigdsq, DblArray *rsq);
#endif
double      Tabulation_sigdsq(Tabulation * tab, int bin, double *sep_cm,
                              long unsigned *nobs);
double      Tabulation_rsq(Tabulation * tab, int bin, double *sep_cm,
                           long unsigned *nobs);
double      Tabulation_denom(Tabulation * tab, int bin, double *sep_cm,
                             long unsigned *nobs);
long unsigned Tabulation_rawCounts(Tabulation * tab, int bin,
                                   double *numerator, double *denominator,
                                   double *sumRsq, double *sep_cm);
long unsigned Tabulation_nObs(const Tabulation * tab, int i);
int         Tabulation_equals(const Tabulation * x, const Tabulation * y);
void        Tabulation_dump(const Tabulation * tab, FILE * ofp);
Tabulation *Tabulation_restore(FILE * ifp);
int         Tabulation_isfinite(const Tabulation * tab);
void        Tabulation_sanityCheck(Tabulation * tab, const char *file,
                                   int line);
int         Tabulation_overflow(Tabulation * tab);

#ifndef NDEBUG
void        Tabulation_test(int verbose);
#endif

/* inline functions */
static inline void Tabulation_record(Tabulation * tab, double Dsq,
                                     double pqpq, double sep_cm,
                                     unsigned wgt);
static inline int getBin(double sep_cm, double bins_per_cm);

static inline int getBin(double sep_cm, double bins_per_cm) {
    return ((int) floor(sep_cm * bins_per_cm));
}

/*
 * Record a single comparison in the appropriate entry of the arrays.
 * Wgt is the number of copies of this comparison. May be 0 or >1 in
 * bootstrap samples.
 */
static inline void Tabulation_record(Tabulation * tab, double Dsq,
                                     double pqpq, double sep_cm,
                                     unsigned wgt) {

    myassert(Dsq >= 0.0);
    myassert(pqpq > 0.0);

    /* determine bin from sep_cm */
    register int ndx = getBin(sep_cm, tab->bins_per_cm);

#ifndef NDEBUG
    if(ndx < 0 || ndx >= tab->nbins) {
        fprintf(stderr, "ERR @%s:%d: ndx=%d; not in [0,%d].\n",
                __FILE__, __LINE__, ndx, tab->nbins - 1);

        if(sep_cm < 0 || sep_cm >= tab->windowcm)
            fprintf(stderr, "Cause: sep_cm=%1.15lg is not in [0,%1.15lg)?\n",
                    sep_cm, tab->windowcm);
        else
            fprintf(stderr, "However, sep_cm=%1.15lg IS in [0,%1.15lg)?\n",
                    sep_cm, tab->windowcm);
        exit(EXIT_FAILURE);
    }
#endif

    tab->nobs[ndx] += wgt;
    tab->numerator[ndx] += wgt * Dsq;
    tab->denominator[ndx] += wgt * pqpq;
    tab->rsq[ndx] += wgt * Dsq / pqpq;
    tab->sep_cm[ndx] += wgt * sep_cm;
}

#endif
