/**
 * @file spectab.c
 * @author Alan R. Rogers
 * @brief Class Spectab tabulates data used to calculate the site frequency spectrum. 
 * @copyright Copyright (c) 2015, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "spectab.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>

struct Spectab {
    unsigned    nSamp;          // haploid sample size
    int         folded;         // 1: folded spectrum; 0: unfolded
    long unsigned nobs;         // number of observations tabulated

    // Unfolded spectrum: s[i] is number of SNPs at which i+1 copies
    // of the derived allele are present in the sample. i ranges from
    // 0 to nSamp-2, so dim is nSamp-1.
    //
    // Folded spectrum: s[i] is number of SNPs at which i+1 copies
    // of the minor allele are present in the sample. i ranges from
    // 0 to floor(nSamp/2) - 1, so dim is floor(nSamp/2).

    unsigned    dim;            // dimension of array s
    long unsigned *s;
};

/**
 * Allocate a new object of type Spectab.
 *
 * @param[in] nSamp The haploid sample size--twice the number of diploid individuals.
 *
 * @param[in] folded 1 for a folded spectrum, 0 for unfolded.
 *
 * @returns A newly allocated object of type Spectab.
 */
Spectab    *Spectab_new(unsigned nSamp, int folded) {
    Spectab *tab = (Spectab *) malloc(sizeof(Spectab));
    checkmem(tab, __FILE__, __LINE__);
    memset(tab, 0, sizeof(Spectab));

    tab->nSamp = nSamp;
    tab->folded = folded;
    switch(folded) {
    case 0:
        tab->dim = nSamp-1;
        break;
    case 1:
        tab->dim = nSamp/2u;
        break;
    default:
        eprintf("%s:%s:%d: folded=%u; must be 0 or 1",
                __FILE__, __func__, __LINE__, folded);
    }

    tab->s = malloc(tab->dim * sizeof(tab->s[0]));
    checkmem(tab->s, __FILE__, __LINE__);
    memset(tab->s, 0, tab->dim * sizeof(tab->s[0]));

#ifndef NDEBUG
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
#endif

    return tab;
}

/*
 * Allocate and initialize a duplicate of a Spectab structure.
 */
Spectab *Spectab_dup(Spectab * old) {
    assert(old);
    Spectab *new = memdup(old, sizeof(Spectab));
    checkmem(new, __FILE__, __LINE__);

    new->s = memdup(old->s, old->dim * sizeof(old->s[0]));
    checkmem(new, __FILE__, __LINE__);

#ifndef NDEBUG
    Spectab_sanityCheck(new, __FILE__, __LINE__);
#endif

    return new;
}

void Spectab_sanityCheck(Spectab * tab, const char *file, int line) {
    REQUIRE(tab != NULL, file, line);
    REQUIRE(tab->nSamp < 10000, file, line);
    REQUIRE(tab->nobs < ULONG_MAX, file, line);
    if(tab->folded)
        REQUIRE(tab->dim == tab->nSamp/2u, file, line);
    else
        REQUIRE(tab->dim + 1 == tab->nSamp, file, line);

    unsigned i;
    unsigned long tot=0;
    for(i=0; i < tab->dim; ++i)
        tot += tab->s[i];
    REQUIRE(tot == tab->nobs, file, line);
    return;
}

/* Return 1 if the two Spectab structs are idential; zero otherwise */
int Spectab_equals(const Spectab * x, const Spectab * y) {

    if(x == NULL || y == NULL)
        return 0;

    if(x->nSamp != y->nSamp
       || x->folded != y->folded
       || x->nobs != y->nobs
       || x->dim != y->dim)
        return 0;

    if(memcmp(x->s, y->s, x->dim * sizeof(x->s[0])))
        return 0;

    return 1;
}

void Spectab_print(Spectab * tab, FILE * ofp) {
    assert(tab);
    fprintf(ofp,
            "Spectab: nSamp  = %u\n"
            "         folded = %d\n"
            "         nobs   = %lu\n"
            "         dim    = %u\n",
            tab->nSamp, tab->folded, tab->nobs, tab->dim);
    fprintf(ofp, "%10s %10s\n", "i", "spec[i]");
    for(unsigned i = 0; i < tab->dim; ++i)
        fprintf(ofp, "%10u %10lu\n", i, tab->s[i]);
}

void Spectab_free(Spectab * tab) {
    if(tab == NULL)
        eprintf("%s:%s:%d: can't free a NULL pointer",
                __FILE__, __func__, __LINE__);
    free(tab);
    free(tab->s);
    return;
}

/// Get data from Spectab object. Function returns nobs.
long unsigned Spectab_report(const Spectab * tab, unsigned dim, int *folded,
                             long unsigned *spec) {
    assert(tab);
    assert(folded);
    assert(spec);
    if(dim != tab->dim)
        eprintf("%s:%s:%d: dim=%u doesn't match tab->dim=%u",
                __FILE__, __func__, __LINE__, dim, tab->dim);

    *folded = tab->folded;

    for(unsigned i = 0; i < tab->dim; ++i)
        spec[i] = tab->s[i];

    return tab->nobs;
}

/*
 * Add Spectab structure y to Spectab structure x and return a pointer
 * to x. It is an error for x to be NULL. If y is NULL, the function
 * is a NOOP.
 */
void Spectab_plus_equals(Spectab * x, const Spectab * y) {
    if(x == NULL) 
        eprintf("%s:s:d: destination is NULL\n",
                __FILE__,__func__,__LINE__);
    if(y == NULL)
        return;
    if(x->dim != y->dim || x->folded != y->folded) 
        eprintf("%s:s:d: unconformable arguments\n",
                __FILE__,__func__,__LINE__);
    if(x == y) 
        eprintf("%s:s:d: x and y are same pointer\n",
                __FILE__,__func__,__LINE__);

    x->nobs += y->nobs;
    for(unsigned i = 0; i < x->dim; ++i)
        x->s[i] += y->s[i];
    return;
}

void Spectab_dump(const Spectab * tab, FILE * ofp) {

    if(0 > fprintf(ofp," %u %d %lu %u\n",
                   tab->nSamp, tab->folded, tab->nobs, tab->dim))
        eprintf("%s:%s:%d: fprintf returned %d instead of 4",
                __FILE__, __func__, __LINE__);

    for(unsigned i = 0; i < tab->dim; ++i) {
        if(0 > fprintf(ofp, " %lu", tab->s[i]))
            eprintf("%s:%s:%d: fprintf", __FILE__, __func__, __LINE__);
    }
    putc('\n', ofp);
}

Spectab *Spectab_restore(FILE * ifp) {
    unsigned nSamp, dim;
    int      folded;
    unsigned long nobs;

    if(4 != fscanf(ifp, "%u %d %lu %u", &nSamp, &folded, &nobs, &dim))
        eprintf("%s:%s:%d: fscanf", __FILE__, __func__, __LINE__);

    Spectab *tab = Spectab_new(nSamp, folded);
    checkmem(tab, __FILE__, __LINE__);
    tab->nobs = nobs;
    if(tab->dim != dim)
        eprintf("%s:%s:%d: tab->dim=%u should match dim=%u",
                __FILE__, __func__, __LINE__, tab->dim, dim);

    for(unsigned i = 0; i < tab->dim; ++i) {
        if(1 != fscanf(ifp, "%lu", tab->s + i))
            eprintf("%s:%s:%d: fscanf", __FILE__, __func__, __LINE__);
    }

#ifndef NDEBUG
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
#endif

    return tab;
}

/*
 * Record a single SNP.  Wgt is the number of copies of this
 * SNP. May be 0 or >1 in bootstrap samples.
 */
void Spectab_record(Spectab * tab, unsigned alleleCount,
                                  unsigned wgt) {
    if(alleleCount == 0 || alleleCount > tab->dim)
        eprintf("%s:%s:%d: Bad alleleCount:%u."
                " Must be in [%u, %u].\n",
                __FILE__, __func__, __LINE__, alleleCount, 0u, tab->dim);

    tab->s[alleleCount - 1] += wgt;
    tab->nobs += wgt;
}



/// Return the number of observations.
long unsigned Spectab_nObs(const Spectab * tab) {
    assert(tab);

    return tab->nobs;
}

#ifndef NDEBUG
void Spectab_test(int verbose) {
    unsigned nSamp = 60;
    int folded;

    /* test Spectab_new and Spectab_free */
    folded = 0;
    Spectab *tab = Spectab_new(nSamp, folded);
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
    Spectab_free(tab);

    folded = 1;
    tab = Spectab_new(nSamp, folded);
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
    Spectab_free(tab);
    unitTstResult("Spectab_new", "OK");
    unitTstResult("Spectab_free", "OK");

    tab = Spectab_new(nSamp, folded);
    Spectab *tab2 = Spectab_dup(tab);
    Spectab_sanityCheck(tab2, __FILE__, __LINE__);
    assert(Spectab_equals(tab, tab2));
    unitTstResult("Spectab_dup", "OK");

    Spectab_record(tab2, nSamp/3, 1u);
    assert(!Spectab_equals(tab, tab2));
    Spectab_free(tab2);
    tab2 = Spectab_new(nSamp, !folded);
    assert(!Spectab_equals(tab, tab2));
    Spectab_free(tab2);
    tab2 = Spectab_new(nSamp-1, folded);
    assert(!Spectab_equals(tab, tab2));
    unitTstResult("Spectab_equals", "OK");

    Spectab_free(tab2);
    tab2 = Spectab_new(nSamp, folded);
    Spectab_record(tab2, 1, 1u); // nobs=1
    Spectab_record(tab2, 1, 1u); // nobs=2
    Spectab_record(tab2, 2, 2u); // nobs=4

    unsigned dim = (folded ? nSamp/2 : nSamp-1);
    int folded2;
    long unsigned nobs;
    long unsigned spec[dim];
    nobs = Spectab_report(tab2, dim, &folded2, spec);
    assert(folded == folded2);
    assert(nobs = 4);
    assert(nobs = Spectab_nObs(tab2));
    assert(spec[0] == 2);
    assert(spec[1] == 2);
    assert(spec[2] == 0);
    unitTstResult("Spectab_report", "OK");
    unitTstResult("Spectab_newSNP", "OK");

    /* Spectab_plus_equals */
    Spectab_plus_equals(tab, tab2);
    assert(tab != tab2);
    assert(Spectab_equals(tab, tab2));
    unitTstResult("Spectab_plus_equals", "OK");

    Spectab_free(tab2);
    tab2 = Spectab_new(nSamp, folded);
    
    const char *fname = "spectab-dummy.tab";
    FILE       *fp = fopen(fname, "w");

    Spectab_dump(tab, fp);
    fclose(fp);
    fp = fopen(fname, "r");
    tab2 = Spectab_restore(fp);
    fclose(fp);
    assert(Spectab_equals(tab, tab2));
    remove(fname);

    unitTstResult("Spectab_dump", "OK");
    unitTstResult("Spectab_restore", "OK");

    if(verbose) {
        Spectab_print(tab, stdout);
        unitTstResult("Spectab_print", "OK");
    }

    Spectab_free(tab);
    Spectab_free(tab2);

    unitTstResult("Spectab", "OK");
}
#endif
