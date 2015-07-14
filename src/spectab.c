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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>

struct Spectab {
    unsigned    twoNsmp;       // haploid sample size
    int         folded;         // true: folded spectrum; false: unfolded
    long unsigned nobs;         // number of observations tabulated

    // Unfolded spectrum: s[i] is number of SNPs at which i+1 copies
    // of the derived allele are present in the sample. i ranges from
    // 0 to twoNsmp-2, so dim is twoNsmp-1.
    //
    // Folded spectrum: s[i] is number of SNPs at which i+1 copies
    // of the minor allele are present in the sample. i ranges from
    // 0 to floor(twoNsmp/2) - 1, so dim is floor(twoNsmp/2).

    unsigned    dim;            // dimension of array s
    long unsigned *s;
};

/// Return dimension of spectrum, based on haploid sample size and
/// the boolean variable "folded", which indicates whether we're
/// dealing with a full site-frequency spectrum or a folded one.
unsigned specdim(unsigned twoNsmp, int folded) {
    return (folded ? twoNsmp/2u : twoNsmp-1u);
}

/**
 * Allocate a new object of type Spectab.
 *
 * @param[in] twoNsmp The haploid sample size--twice the number of
 * diploid individuals.  
 *
 * @param[in] folded 1 for a folded spectrum, 0 for unfolded.
 *
 * @returns A newly allocated object of type Spectab.
 */
Spectab    *Spectab_new(unsigned twoNsmp, int folded) {
    Spectab *tab = (Spectab *) malloc(sizeof(Spectab));
    CHECKMEM(tab);
    memset(tab, 0, sizeof(Spectab));

    tab->twoNsmp = twoNsmp;
    tab->folded = folded;
    tab->dim = specdim(twoNsmp, folded);
    tab->s = malloc(tab->dim * sizeof(tab->s[0]));
    CHECKMEM(tab->s);
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
    CHECKMEM(new);

    new->s = memdup(old->s, old->dim * sizeof(old->s[0]));
    CHECKMEM(new);

#ifndef NDEBUG
    Spectab_sanityCheck(new, __FILE__, __LINE__);
#endif

    return new;
}

void Spectab_sanityCheck(Spectab * tab, const char *file, int line) {
    REQUIRE(tab != NULL, file, line);
    REQUIRE(tab->twoNsmp < 10000, file, line);
    REQUIRE(tab->nobs < ULONG_MAX, file, line);
    if(tab->folded)
        REQUIRE(tab->dim == tab->twoNsmp/2u, file, line);
    else
        REQUIRE(tab->dim + 1 == tab->twoNsmp, file, line);

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
        return false;

    if(x->twoNsmp != y->twoNsmp
       || x->folded != y->folded
       || x->nobs != y->nobs
       || x->dim != y->dim)
        return false;

    if(memcmp(x->s, y->s, x->dim * sizeof(x->s[0])))
        return false;

    return true;
}

void Spectab_print(Spectab * tab, FILE * ofp) {
    assert(tab);
    fprintf(ofp,
            "Spectab: twoNsmp  = %u\n"
            "         folded    = %d\n"
            "         nobs      = %lu\n"
            "         dim       = %u\n",
            tab->twoNsmp, tab->folded, tab->nobs, tab->dim);
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

int         Spectab_folded(Spectab *self) {
    return self->folded;
}

/// Get spectrum from Spectab object. Function returns nobs. 
long unsigned Spectab_report(const Spectab * tab, ULIntArray *spec) {
    assert(tab);
    assert(spec);
#ifndef NDEBUG
    if(tab->dim != ULIntArray_dim(spec))
        eprintf("%s:%s:%d: dimemsion mismatch", __FILE__, __func__, __LINE__);
#endif

    for(unsigned i = 0; i < tab->dim; ++i)
        ULIntArray_set(spec, i, tab->s[i]);

    return tab->nobs;
}

/// Return the number of sites at which the counted allele is present i+1 times
long unsigned Spectab_get(const Spectab *st, unsigned i) {
    assert(i >= 0);
    assert(i < st->dim);
    return st->s[i];
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
                   tab->twoNsmp, tab->folded, tab->nobs, tab->dim))
        eprintf("%s:%s:%d: fprintf returned %d instead of 4",
                __FILE__, __func__, __LINE__);

    for(unsigned i = 0; i < tab->dim; ++i) {
        if(0 > fprintf(ofp, " %lu", tab->s[i]))
            eprintf("%s:%s:%d: fprintf", __FILE__, __func__, __LINE__);
    }
    putc('\n', ofp);
}

Spectab *Spectab_restore(FILE * ifp) {
    unsigned twoNsmp, dim;
    int      folded;
    unsigned long nobs;

    if(4 != fscanf(ifp, "%u %d %lu %u", &twoNsmp, &folded, &nobs, &dim))
        eprintf("%s:%s:%d: fscanf", __FILE__, __func__, __LINE__);

    Spectab *tab = Spectab_new(twoNsmp, folded);
    CHECKMEM(tab);
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
#ifndef NDEBUG
    if(alleleCount == 0 || alleleCount > tab->dim)
        eprintf("%s:%s:%d: Bad alleleCount:%u."
                " Must be in [%u, %u].\n",
                __FILE__, __func__, __LINE__, alleleCount, 0u, tab->dim);
#endif

    tab->s[alleleCount - 1] += wgt;
    tab->nobs += wgt;
}



/// Return the number of observations.
long unsigned Spectab_nObs(const Spectab * tab) {
    assert(tab);

    return tab->nobs;
}

/// Return dimension of internal array
unsigned    Spectab_dim(const Spectab * st) {
    return st->dim;
}


#ifndef NDEBUG
void Spectab_test(int verbose) {
    unsigned twoNsmp = 60;
    int folded;

    REQUIRE(twoNsmp/2u == specdim(twoNsmp, true), __FILE__, __LINE__);
    REQUIRE(twoNsmp-1 == specdim(twoNsmp, false), __FILE__, __LINE__);
    unitTstResult("specdim", "OK");

    /* test Spectab_new and Spectab_free */
    folded = false;
    unsigned dim = (folded ? twoNsmp/2 : twoNsmp-1);
    Spectab *tab = Spectab_new(twoNsmp, folded);
    assert(dim == Spectab_dim(tab));
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
    Spectab_free(tab);

    folded = true;
    tab = Spectab_new(twoNsmp, folded);
    Spectab_sanityCheck(tab, __FILE__, __LINE__);
    Spectab_free(tab);
    unitTstResult("Spectab_new", "OK");
    unitTstResult("Spectab_free", "OK");

    tab = Spectab_new(twoNsmp, folded);
    Spectab *tab2 = Spectab_dup(tab);
    Spectab_sanityCheck(tab2, __FILE__, __LINE__);
    assert(Spectab_equals(tab, tab2));
    unitTstResult("Spectab_dup", "OK");

    Spectab_record(tab2, twoNsmp/3, 1u);
    assert(!Spectab_equals(tab, tab2));
    Spectab_free(tab2);
    tab2 = Spectab_new(twoNsmp, !folded);
    assert(!Spectab_equals(tab, tab2));
    Spectab_free(tab2);
    tab2 = Spectab_new(twoNsmp-1, folded);
    assert(!Spectab_equals(tab, tab2));
    unitTstResult("Spectab_equals", "OK");

    Spectab_free(tab2);
    tab2 = Spectab_new(twoNsmp, folded);
    Spectab_record(tab2, 1, 1u); // nobs=1
    Spectab_record(tab2, 1, 1u); // nobs=2
    Spectab_record(tab2, 2, 2u); // nobs=4

    dim = (folded ? twoNsmp/2 : twoNsmp-1);
    long unsigned nobs;
    ULIntArray *spec = ULIntArray_new(dim);
    CHECKMEM(spec);

    nobs = Spectab_report(tab2, spec);
    assert(dim == Spectab_dim(tab2));
    assert(nobs = 4);
    assert(nobs = Spectab_nObs(tab2));
    assert(2 == ULIntArray_get(spec, 0));
    assert(2 == ULIntArray_get(spec, 1));
    assert(0 == ULIntArray_get(spec, 2));
    assert(2 == Spectab_get(tab2, 0));
    assert(2 == Spectab_get(tab2, 1));
    assert(0 == Spectab_get(tab2, 2));
    unitTstResult("Spectab_report", "OK");
    unitTstResult("Spectab_newSNP", "OK");

    /* Spectab_plus_equals */
    Spectab_plus_equals(tab, tab2);
    assert(tab != tab2);
    assert(Spectab_equals(tab, tab2));
    unitTstResult("Spectab_plus_equals", "OK");

    Spectab_free(tab2);
    tab2 = Spectab_new(twoNsmp, folded);
    
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

    ULIntArray_free(spec);
    Spectab_free(tab);
    Spectab_free(tab2);

    unitTstResult("Spectab", "OK");
}
#endif
