/**
 * @file tabulation.c
 * @author Alan R. Rogers
 * @brief Class Tabulation tabulates data used to calculate \f$sigma_d^2\f$.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "tabulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

static void Tabulation_initArrays(Tabulation * tab);
static void tstGetBin(double sep_cm, double windowcm, int trubin, int nbins,
                      const char *file, int line);
static size_t Tabulation_size(int nbins);
static inline double get_bins_per_cm(int nbins, double windowcm);

static void tstGetBin(double sep_cm, double windowcm, int trubin, int nbins,
                      const char *file, int line) {

    double      bins_per_cm = get_bins_per_cm(nbins, windowcm);
    int         bin = getBin(sep_cm, bins_per_cm);

    if(bin != trubin)
        eprintf("ERR@%s:%d: getBin(%lf,%lf)=%d rather than %d",
                file, line, sep_cm, bins_per_cm, bin, trubin);
}

static inline double get_bins_per_cm(int nbins, double windowcm) {
    return nbins / windowcm;
}

size_t Tabulation_size(int nbins) {
    Tabulation  dummy;
    size_t      size = sizeof(Tabulation);

    size += nbins * sizeof(dummy.nobs[0]);
    size += nbins * sizeof(dummy.sep_cm[0]);
    size += nbins * sizeof(dummy.numerator[0]);
    size += nbins * sizeof(dummy.denominator[0]);
    size += nbins * sizeof(dummy.rsq[0]);

    return size;
}

static void Tabulation_initArrays(Tabulation * tab) {
    assert(tab);

    tab->nobs = (void *) &(tab[1]);
    tab->sep_cm = (void *) &(tab->nobs[tab->nbins]);
    tab->numerator = (void *) &(tab->sep_cm[tab->nbins]);
    tab->denominator = (void *) &(tab->numerator[tab->nbins]);
    tab->rsq = (void *) &(tab->denominator[tab->nbins]);

#ifndef NDEBUG
    /*
     * tst: the size of the chunk spanned by all arrays should equal size.
     */
    size_t      tst = ((size_t) & (tab->rsq[tab->nbins])
                       - ((size_t) tab));

    if(Tabulation_size(tab->nbins) != tst) {
        fprintf(stderr, "\nERR@%s:%d:"
                " size=%lu != tst=%lu\n\n",
                __FILE__, __LINE__,
                (long unsigned) Tabulation_size(tab->nbins),
                (long unsigned) tst);
        dostacktrace(__FILE__, __LINE__, stderr);
        exit(EXIT_FAILURE);
    }
#endif
}

/**
 * Allocate a new object of type Tabulation.
 *
 * This is accomplished via a single call to malloc, so that the
 * entire object will reside in contiguous memory. This minimizes page
 * faults and should make time-critical functions such as
 * Tabulation_record run faster.
 *
 * @param[in] windowcm The size in centimorgans of the window used to scan
 * the genome.
 *
 * @param[in] nbins The number of bins in which to tabulate values.
 *
 * @returns A newly allocated object of type Tabulation, which should
 * be freed with Tabulation_free.
 */
Tabulation *Tabulation_new(double windowcm, int nbins) {
    size_t      size = Tabulation_size(nbins);
    Tabulation *tab = (Tabulation *) malloc(size);

    checkmem(tab, __FILE__, __LINE__);
    memset(tab, 0, size);

    tab->nbins = nbins;
    tab->windowcm = windowcm;
    tab->bins_per_cm = nbins / windowcm;

    Tabulation_initArrays(tab);

#ifndef NDEBUG
    Tabulation_sanityCheck(tab, __FILE__, __LINE__);
#endif

    return tab;
}

/*
 * Allocate and initialize a duplicate of a Tabulation structure. May
 * be freed with Tabulation_free.
 */
Tabulation *Tabulation_dup(Tabulation * old) {
    assert(old);
    Tabulation *new = memdup(old, Tabulation_size(old->nbins));

    Tabulation_initArrays(new);

#ifndef NDEBUG
    Tabulation_sanityCheck(new, __FILE__, __LINE__);
#endif

    return new;
}

void Tabulation_sanityCheck(Tabulation * tab, const char *file, int line) {
    if(tab == NULL) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: NULL pointer", file, __func__, line);
    }
    if(tab->nbins <= 0) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: nbins=%d", file, __func__, line, tab->nbins);
    }
    if(tab->windowcm <= 0 || tab->windowcm > 1.0) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: windowcm=%lg", file, __func__, line,
                tab->windowcm);
    }
    int         nbins2 = 0.5 + tab->bins_per_cm * tab->windowcm;

    if(nbins2 != tab->nbins) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: nbins2=%d != nbins=%d", file, __func__, line,
                nbins2, tab->nbins);
    }

    if(tab->nobs != (void *) ((size_t) tab + sizeof(Tabulation))) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: nobs incorrect", file, __func__, line);
    }
    if(tab->sep_cm != (void *) ((size_t) tab->nobs
                                + tab->nbins * sizeof(long))) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: sep_cm incorrect", file, __func__, line);
    }
    if(tab->numerator != (void *) ((size_t) tab->sep_cm
                                   + tab->nbins * sizeof(double))) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: numerator incorrect", file, __func__, line);
    }
    if(tab->denominator != (void *) ((size_t) tab->numerator
                                     + tab->nbins * sizeof(double))) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: denominator incorrect", file, __func__, line);
    }
    if(tab->rsq != (void *) ((size_t) tab->denominator
                             + tab->nbins * sizeof(double))) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: rsq incorrect", file, __func__, line);
    }

    void       *end = (void *) ((size_t) tab + Tabulation_size(tab->nbins));

    if((void *) &tab->nobs[tab->nbins] > end) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: nobs end incorrect", file, __func__, line);
    }
    if((void *) &tab->sep_cm[tab->nbins] > end) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: sep_cm end incorrect", file, __func__, line);
    }
    if((void *) &tab->numerator[tab->nbins] > end) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: numerator end incorrect", file, __func__,
                line);
    }
    if((void *) &tab->denominator[tab->nbins] > end) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: denominator end incorrect", file, __func__,
                line);
    }
    if((void *) &tab->rsq[tab->nbins] > end) {
        dostacktrace(file, line, stderr);
        eprintf("ERR:%s:%s:%d: rsq end incorrect", file, __func__, line);
    }
    return;
}

/** Check for overflow */
int Tabulation_overflow(Tabulation * tab) {
    if(tab->nobs[0] < tab->oldNobs[0])
        return 1;
    tab->oldNobs[0] = tab->nobs[0];

    if(tab->nobs[tab->nbins - 1] < tab->oldNobs[1])
        return 1;
    tab->oldNobs[1] = tab->nobs[tab->nbins - 1];

    return 0;
}

/* Return 1 if the two Tabulation structs are idential; zero otherwise */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
int Tabulation_equals(const Tabulation * x, const Tabulation * y) {

    if(x == NULL || y == NULL)
        return 0;

    if(x->nbins != y->nbins || x->windowcm != y->windowcm)
        return 0;

    if(memcmp(x->nobs, y->nobs, x->nbins * sizeof(x->nobs[0])))
        return 0;

    if(memcmp(x->sep_cm, y->sep_cm, x->nbins * sizeof(x->sep_cm[0])))
        return 0;

    if(memcmp(x->numerator, y->numerator, x->nbins * sizeof(x->numerator[0])))
        return 0;

    if(memcmp(x->denominator, y->denominator,
              x->nbins * sizeof(x->denominator[0])))
        return 0;

    if(memcmp(x->rsq, y->rsq, x->nbins * sizeof(x->rsq[0])))
        return 0;

    return 1;
}
#pragma GCC diagnostic pop

void Tabulation_print(Tabulation * tab, FILE * ofp) {
    int         i;

    myassert(tab);

    fprintf(ofp, "Tabulation: windowcm = %lg\n", tab->windowcm);
    fprintf(ofp, "%10s %10s %10s %10s %10s\n",
            "nobs", "sep", "num", "denom", "sumRsq");
    for(i = 0; i < tab->nbins; ++i)
        fprintf(ofp, "%10lu %10g %10g %10g %10g\n", tab->nobs[i],
                tab->sep_cm[i], tab->numerator[i], tab->denominator[i],
                tab->rsq[i]);
}

void Tabulation_free(Tabulation * tab) {
    if(tab == NULL)
        return;
    free(tab);
    return;
}

/**
 * Return 1 if tab is free of non-finite values, 0 otherwise.
 */
int Tabulation_isfinite(const Tabulation * tab) {
    int         i;

    myassert(tab);

    if(!isfinite(tab->windowcm + tab->bins_per_cm))
        return 0;

    for(i = 0; i < tab->nbins; ++i) {
        if(!isfinite(tab->sep_cm[i]
                     + tab->numerator[i]
                     + tab->denominator[i]
                     + tab->rsq[i]))
            return 0;
    }
    return 1;
}

/*
 * Fill arrays sep_cm, nobs, sigdsq, and rsq with estimates from the
 * tabulation. It is legal for nobs and/or rsq to equal NULL, in which
 * case no values are stored. Returns 0 on success, or 1 if data are invalid
 * because an overflow has occurred.
 */
#if 0
int Tabulation_report(Tabulation * tab, double *sep_cm, long unsigned *nobs,
                              double *sigdsq, double *rsq) {
    int         i;
    unsigned    n;

    myassert(tab);
    myassert(sep_cm);
    myassert(sigdsq);

    if(Tabulation_overflow(tab))
        return 1;

    for(i = 0; i < tab->nbins; ++i) {
        n = tab->nobs[i];
        if(nobs)
            nobs[i] = n;

        if(n == 0) {
            sep_cm[i] = sigdsq[i] = strtod("NAN", NULL);
            if(rsq)
                rsq[i] = strtod("NAN", NULL);
        } else {
            sigdsq[i] = tab->numerator[i] / tab->denominator[i];
            sep_cm[i] = tab->sep_cm[i] / n;
            if(rsq)
                rsq[i] = tab->rsq[i] / n;
        }
    }

    return 0;
}
#else
int         Tabulation_report(Tabulation *tab,
                              DblArray *sep_cm, ULIntArray *nobs,
                              DblArray *sigdsq, DblArray *rsq) {
    int         i;
    unsigned long  n;

    myassert(tab);
    myassert(sep_cm);
    myassert(sigdsq);
    assert(tab->nbins == DblArray_dim(sep_cm));
    assert(tab->nbins == DblArray_dim(sigdsq));
    assert(tab->nbins == DblArray_dim(rsq));
    assert(tab->nbins == ULIntArray_dim(nobs));

    if(Tabulation_overflow(tab))
        return 1;

    for(i = 0; i < tab->nbins; ++i) {
        n = tab->nobs[i];
        if(nobs)
            ULIntArray_set(nobs, i, n);

        if(n == 0) {
            double nanval = strtod("NAN", NULL);
            DblArray_set(sep_cm, i, nanval);
            DblArray_set(sigdsq, i, nanval);
            if(rsq)
                DblArray_set(rsq, i, nanval);
        } else {
            DblArray_set(sigdsq, i, tab->numerator[i] / tab->denominator[i]);
            DblArray_set(sep_cm, i, tab->sep_cm[i] / n);
            if(rsq)
                DblArray_set(rsq, i, tab->rsq[i] / n);
        }
    }

    return 0;
}
#endif

/*
 * Return the sigdsq estimate for a single bin. On entry, "tab"
 * provides the data to be used, and "bin" specifies the bin in which
 * we are interested. On return, *sep_cm will equal the mean
 * separation between the pairs of sites tabulated, *nobs will equal
 * the number of observations in that bin, and the function will
 * return the value of sigma_d^2 estimated from this bin. If there are
 * no observations in the bin, the function returns NAN.
 */
double Tabulation_sigdsq(Tabulation * tab, int bin, double *sep_cm,
                         long unsigned *nobs) {
    double      sigdsq;

    myassert(tab);
    myassert(sep_cm);
    myassert(nobs);
    myassert(bin < tab->nbins);
    myassert(bin >= 0);
    *nobs = tab->nobs[bin];
    if(tab->nobs[bin] == 0)
        return strtod("NAN", NULL);

    sigdsq = tab->numerator[bin] / tab->denominator[bin];
    *sep_cm = tab->sep_cm[bin] / tab->nobs[bin];

    return sigdsq;
}

/*
 * Return the rsq estimate for a single bin. On entry, "tab"
 * provides the data to be used, and "bin" specifies the bin in which
 * we are interested. On return, *sep_cm will equal the mean
 * separation between the pairs of sites tabulated, *nobs will equal
 * the number of observations in that bin, and the function will
 * return the value of rsq estimated from this bin. If there are
 * no observations in the bin, the function returns NAN.
 */
double Tabulation_rsq(Tabulation * tab, int bin, double *sep_cm,
                      long unsigned *nobs) {
    double      rsq;

    myassert(tab);
    myassert(sep_cm);
    myassert(nobs);
    myassert(bin < tab->nbins);
    myassert(bin >= 0);
    *nobs = tab->nobs[bin];
    if(tab->nobs[bin] == 0)
        return strtod("NAN", NULL);

    rsq = tab->rsq[bin] / tab->nobs[bin];
    *sep_cm = tab->sep_cm[bin] / tab->nobs[bin];

    return rsq;
}

/*
 * Return the raw counts for a single bin. Counts for numerator,
 * denominator, sep_cm, and sumRsq are returned in the corresponding
 * arguments. The function itself returns nobs.
 */
long unsigned Tabulation_rawCounts(Tabulation * tab, int bin,
                                   double *numerator, double *denominator,
                                   double *sumRsq, double *sep_cm) {
    myassert(tab);
    myassert(numerator);
    myassert(denominator);
    myassert(sep_cm);
    myassert(bin < tab->nbins);
    myassert(bin >= 0);

    *numerator = tab->numerator[bin];
    *denominator = tab->denominator[bin];
    *sep_cm = tab->sep_cm[bin];
    *sumRsq = tab->rsq[bin];

    return tab->nobs[bin];
}

/**
 * Return the number of observations in the i'th bin.
 */
long unsigned Tabulation_nObs(const Tabulation * tab, int i) {
    myassert(tab);
    myassert(i < tab->nbins);

    return tab->nobs[i];
}

/*
 * Add Tabulation structure y to Tabulation structure x and return a pointer
 * to x. It is an error for x to be NULL. If y is NULL, the function
 * is a NOOP.
 */
void Tabulation_plus_equals(Tabulation * x, const Tabulation * y) {
    int         i;

    if(x == NULL) {
        fprintf(stderr, "%s: destination is NULL\n", __func__);
        exit(1);
    }
    if(y == NULL) {
        fprintf(stderr, "%s: source is NULL\n", __func__);
        return;
    }
    if(x->nbins != y->nbins || x->windowcm != y->windowcm) {
        fprintf(stderr, "%s: unconformable arguments\n", __func__);
        exit(1);
    }
    if(x == y) {
        fprintf(stderr, "%s: x and y are same pointer\n", __func__);
    }
    for(i = 0; i < x->nbins; ++i) {
        x->nobs[i] += y->nobs[i];
        x->sep_cm[i] += y->sep_cm[i];
        x->numerator[i] += y->numerator[i];
        x->denominator[i] += y->denominator[i];
        x->rsq[i] += y->rsq[i];
    }
    return;
}

void Tabulation_dump(const Tabulation * tab, FILE * ofp) {
    int         bin;
    const static char *fmt = " %0.20lf";

    if(0 > fprintf(ofp, "%d", tab->nbins)
       || 0 > fprintf(ofp, fmt, tab->windowcm)
       || 0 > fprintf(ofp, fmt, tab->bins_per_cm)
        )
        eprintf("%s:%s:%d: fprintf", __FILE__, __func__, __LINE__);
    putc('\n', ofp);

    for(bin = 0; bin < tab->nbins; ++bin) {
        if(0 > fprintf(ofp, " %lu", tab->nobs[bin])
           || 0 > fprintf(ofp, fmt, tab->sep_cm[bin])
           || 0 > fprintf(ofp, fmt, tab->numerator[bin])
           || 0 > fprintf(ofp, fmt, tab->denominator[bin])
           || 0 > fprintf(ofp, fmt, tab->rsq[bin])
            )
            eprintf("%s:%s:%d: fprintf", __FILE__, __func__, __LINE__);
        putc('\n', ofp);
    }
}

Tabulation *Tabulation_restore(FILE * ifp) {
    int         bin, rval, nbins;
    size_t      size;
    double      windowcm, bins_per_cm;

    rval = fscanf(ifp, "%d %lf %lf", &nbins, &windowcm, &bins_per_cm);
    if(rval != 3)
        eprintf("ERR@%s:%d: fscanf returned %d instead of 3",
                __FILE__, __LINE__, rval);

    size = Tabulation_size(nbins);
    Tabulation *tab = malloc(size);

    checkmem(tab, __FILE__, __LINE__);

    memset(tab, 0, size);
    tab->nbins = nbins;
    tab->windowcm = windowcm;
    tab->bins_per_cm = bins_per_cm;
    Tabulation_initArrays(tab);

    for(bin = 0; bin < tab->nbins; ++bin) {
        rval = fscanf(ifp, "%lu", tab->nobs + bin);
        if(rval != 1)
            eprintf("ERR@%s:%d: fscanf returned %d instead of 1",
                    __FILE__, __LINE__, rval);

        rval = fscanf(ifp, "%lf", tab->sep_cm + bin);
        if(rval != 1)
            eprintf("ERR@%s:%d: fscanf returned %d instead of 1",
                    __FILE__, __LINE__, rval);

        rval = fscanf(ifp, "%lf", tab->numerator + bin);
        if(rval != 1)
            eprintf("ERR@%s:%d: fscanf returned %d instead of 1",
                    __FILE__, __LINE__, rval);

        rval = fscanf(ifp, "%lf", tab->denominator + bin);
        if(rval != 1)
            eprintf("ERR@%s:%d: fscanf returned %d instead of 1",
                    __FILE__, __LINE__, rval);

        rval = fscanf(ifp, "%lf", tab->rsq + bin);
        if(rval != 1)
            eprintf("ERR@%s:%d: fscanf returned %d instead of 1",
                    __FILE__, __LINE__, rval);
    }

#ifndef NDEBUG
    Tabulation_sanityCheck(tab, __FILE__, __LINE__);
#endif

    return tab;
}

#ifndef NDEBUG
void Tabulation_test(int verbose) {
    int ok2;

    int         nbins = 10, bin, i;
    double      windowcm = 100.0, Dsq, pqpq, sep_cm, rsq;

    /* lower edge of 0th bin */
    tstGetBin(0.0, windowcm, 0, nbins, __FILE__, __LINE__);

    /* upper edge of last bin */
    tstGetBin(99.999999, windowcm, nbins - 1, nbins, __FILE__, __LINE__);

    /* 1 more than upper edge of last bin */
    tstGetBin(windowcm, windowcm, nbins, nbins, __FILE__, __LINE__);

    /* upper edge of first bin */
    tstGetBin(9.99999999, windowcm, 0, nbins, __FILE__, __LINE__);

    /* lower edge of second bin */
    tstGetBin(10, windowcm, 1, nbins, __FILE__, __LINE__);

    unitTstResult("getBin", "OK");

    assert(Tabulation_size(10) == sizeof(Tabulation)
           + 10 * (sizeof(long) + 4 * sizeof(double)));

    unitTstResult("Tabulation_size", "OK");

    /* test Tabulation_initArrays */
    ok2 = 1;
    nbins = 10;
    Tabulation *tab = malloc(Tabulation_size(nbins));

    memset(tab, 0, Tabulation_size(nbins));
    tab->nbins = nbins;
    Tabulation_initArrays(tab);
    size_t      tst = ((size_t) & (tab->rsq[tab->nbins])
                       - ((size_t) tab));

    if(Tabulation_size(tab->nbins) != tst) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    if(((size_t) tab->nobs - (size_t) tab) != sizeof(Tabulation)) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    if(((size_t) tab->sep_cm - (size_t) tab->nobs)
       != nbins * sizeof(long)) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    if(((size_t) tab->numerator - (size_t) tab->sep_cm)
       != nbins * sizeof(double)) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    if(((size_t) tab->denominator - (size_t) tab->numerator)
       != nbins * sizeof(double)) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    if(((size_t) tab->rsq - (size_t) tab->denominator)
       != nbins * sizeof(double)) {
        printf("%s:%d:%s: FAILING\n",__FILE__,__LINE__,__func__);
        ok2 = 0;
    }
    free(tab);

    if(ok2)
        unitTstResult("Tabulation_initArrays", "OK");
    else {
        unitTstResult("Tabulation_initArrays", "FAIL");
    }

    /* test Tabulation_new and Tabulation_free */
    windowcm = 0.3;
    nbins = 10;
    tab = Tabulation_new(windowcm, nbins);
    assert(Tabulation_isfinite(tab));
    assert(nbins == tab->nbins);
    assert(windowcm == tab->windowcm);
    assert(fabs(nbins / windowcm - tab->bins_per_cm) < 10 * DBL_EPSILON);
    for(i = 0; i < nbins; ++i)
        assert(0 == Tabulation_nObs(tab, i));

    unitTstResult("Tabulation_nObs", "OK");

    /* check organization of arrays */
    assert(((size_t) tab->nobs) == ((size_t) tab) + sizeof(Tabulation));
    assert(((size_t) tab->sep_cm) == ((size_t) tab)
           + sizeof(Tabulation)
           + nbins * sizeof(tab->nobs[0]));
    assert(((size_t) tab->numerator) == ((size_t) tab)
           + sizeof(Tabulation)
           + nbins * sizeof(tab->nobs[0])
           + nbins * sizeof(tab->sep_cm[0]));
    assert(((size_t) tab->denominator) == ((size_t) tab)
           + sizeof(Tabulation)
           + nbins * sizeof(tab->nobs[0])
           + nbins * sizeof(tab->sep_cm[0])
           + nbins * sizeof(tab->numerator[0]));
    assert(((size_t) tab->rsq) == ((size_t) tab)
           + sizeof(Tabulation)
           + nbins * sizeof(tab->nobs[0])
           + nbins * sizeof(tab->sep_cm[0])
           + nbins * sizeof(tab->numerator[0])
           + nbins * sizeof(tab->denominator[0]));
    assert(((size_t) & tab->rsq[nbins] - ((size_t) tab))
           == Tabulation_size(nbins));
    assert((((size_t) tab->rsq)
            + nbins * sizeof(tab->rsq[0])
            - ((size_t) tab)
           ) == Tabulation_size(nbins));

    unitTstResult("Tabulation_new", "OK");

    Tabulation_free(tab);
    unitTstResult("Tabulation_free", "OK");

    /* test Tabulation_record */
    windowcm = 0.3;
    tab = Tabulation_new(windowcm, nbins);

    assert(tab->windowcm == windowcm);
    assert(tab->nbins == nbins);
    assert(fabs(nbins / windowcm - tab->bins_per_cm) < 10 * DBL_EPSILON);

    /* check status of last value in tab. Should be zero */
    size_t      end = ((size_t) tab) + Tabulation_size(nbins);
    double     *xptr = (double *) (end - sizeof(xptr[0]));

    assert(*xptr == 0.0);

    for(i = 0; i < nbins; ++i) {
        assert(tab->nobs[i] == 0);
        assert(tab->sep_cm[i] == 0.0);
        assert(tab->numerator[i] == 0.0);
        if(!((tab->denominator[i] == 0.0)))
            printf("i=%d: ptr=%p; denom[i]=%lg; should be 0\n",
                   i, &tab->denominator[i], tab->denominator[i]);
        assert(tab->denominator[i] == 0.0);
        assert(tab->rsq[i] == 0.0);
    }

    Dsq = 0.25;
    pqpq = 0.75;
    rsq = Dsq / pqpq;
    sep_cm = 0.999999 * windowcm;
    bin = getBin(sep_cm, tab->bins_per_cm);
    assert(bin == nbins - 1);
    Tabulation_record(tab, Dsq, pqpq, sep_cm, 1u);

    for(i = 0; i < nbins; ++i) {
        if(i == bin) {
            assert(tab->nobs[i] == 1);
            assert(tab->sep_cm[i] == sep_cm);
            assert(tab->numerator[i] == Dsq);
            assert(tab->denominator[i] == pqpq);
            assert(tab->rsq[i] == rsq);
        } else {
            assert(tab->nobs[i] == 0);
            assert(tab->sep_cm[i] == 0.0);
            assert(tab->numerator[i] == 0.0);
            assert(tab->denominator[i] == 0.0);
            assert(tab->rsq[i] == 0.0);
        }
    }
    unitTstResult("Tabulation_record", "OK");

    if(verbose)
        Tabulation_print(tab, stdout);

    /* give tab a complicated value  */
    for(i = 0; i < tab->nbins; ++i) {
        tab->nobs[i] = rand();
        tab->numerator[i] = rand();
        tab->denominator[i] = rand();
        tab->sep_cm[i] = rand();
        tab->nobs[i] = rand();
    }

    Tabulation *tab2 = Tabulation_dup(tab);

    assert(tab2 != tab);

    assert(tab2->nobs != tab->nobs);
    assert(tab2->sep_cm != tab->sep_cm);
    assert(tab2->numerator != tab->numerator);
    assert(tab2->denominator != tab->denominator);
    assert(tab2->rsq != tab->rsq);
    assert(Tabulation_equals(tab, tab2));
    unitTstResult("Tabulation_dup", "OK");
    unitTstResult("Tabulation_equals", "OK");

    /* Tabulation_plus_equals */
    Tabulation_plus_equals(tab2, tab);
    assert(tab2 != tab);

    assert(tab->nbins == tab2->nbins);
    assert(tab->windowcm == tab2->windowcm);
    assert(tab->bins_per_cm == tab2->bins_per_cm);
    for(i = 0; i < tab->nbins; ++i) {
        assert(tab2->nobs[i] == 2 * tab->nobs[i]);
        assert(Dbl_near(tab2->sep_cm[i], 2.0 * tab->sep_cm[i]));
        assert(Dbl_near(tab2->numerator[i], 2.0 * tab->numerator[i]));
        assert(Dbl_near(tab2->denominator[i], 2.0 * tab->denominator[i]));
    }
    unitTstResult("Tabulation_plus_equals", "OK");

    /* test Tabulation_sigdsq */
    for(i = 0; i < tab->nbins; ++i) {
        tab->nobs[i] = 10;
        tab->numerator[i] = 1.0;
        tab->denominator[i] = 2.0;
        tab->rsq[i] = 3.0;
        tab->sep_cm[i] = 4.0;
        tab->nobs[i] = 10;
    }

    for(i = 0; i < tab->nbins; ++i) {
        double      sep, num, denom, sumRsq;
        long unsigned n;

        assert(Dbl_near(0.5, Tabulation_sigdsq(tab, i, &sep, &n)));
        assert(0.4 == sep);
        assert(10 == n);
        assert(Dbl_near(0.3, Tabulation_rsq(tab, i, &sep, &n)));
        n = Tabulation_rawCounts(tab, i, &num, &denom, &sumRsq, &sep);
        assert(1.0 == num);
        assert(2.0 == denom);
        assert(3.0 == sumRsq);
        assert(4.0 == sep);
        assert(10 == n);
    }
    unitTstResult("Tabulation_sigdsq", "OK");
    unitTstResult("Tabulation_rawCounts", "OK");

    Tabulation_free(tab);
    tab = NULL;
    Tabulation_free(tab2);
    tab2 = NULL;

    /* give tab a complicated value  */
    windowcm = 0.3;
    nbins = 10;
    tab = Tabulation_new(windowcm, nbins);
    for(i = 0; i < tab->nbins; ++i) {
        tab->nobs[i] = rand();
        tab->numerator[i] = rand();
        tab->denominator[i] = rand();
        tab->rsq[i] = rand();
        tab->sep_cm[i] = rand();
        tab->nobs[i] = rand();
    }

    const char *fname = "dummy.tab";
    FILE       *fp = fopen(fname, "w");

    Tabulation_dump(tab, fp);
    fclose(fp);
    fp = fopen(fname, "r");
    tab2 = Tabulation_restore(fp);
    fclose(fp);
    assert(Tabulation_equals(tab, tab2));
    remove(fname);

    unitTstResult("Tabulation_dump", "OK");
    unitTstResult("Tabulation_restore", "OK");

    Tabulation_free(tab);
    Tabulation_free(tab2);

    unitTstResult("Tabulation", "OK");
}
#endif
