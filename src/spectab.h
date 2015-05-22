/**
 * @file spectab.h
 * @author Alan R. Rogers
 * @brief Header for spectab.c
 * @copyright Copyright (c) 2015, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SPECTAB_H
#  define LDPSIZ_SPECTAB_H

#  include "typedefs.h"
#  include "misc.h"

/*
 * spectab.h. This file defines the Spectab object, which tabulates
 * data for the site frequency spectrum. 
 *
 * The definition of struct Spectab is exposed here (rather than
 * hidden in spectab.c) for the benefit of the inline functions
 * defined at the bottom of this file.
 */
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

Spectab    *Spectab_new(unsigned nSamp, int folded);
Spectab    *Spectab_dup(Spectab * old);
void        Spectab_free(Spectab * tab);
void        Spectab_plus_equals(Spectab * x, const Spectab * y);
void        Spectab_print(Spectab * tab, FILE * ofp);
long unsigned Spectab_report(const Spectab * tab, unsigned dim, int *folded,
                             long unsigned *spec);
int         Spectab_equals(const Spectab * x, const Spectab * y);
void        Spectab_dump(const Spectab * tab, FILE * ofp);
Spectab    *Spectab_restore(FILE * ifp);
void        Spectab_sanityCheck(Spectab * tab, const char *file, int line);
void        Spectab_newSNP(Spectab *tab, unsigned alleleCount);

#  ifndef NDEBUG
void        Spectab_test(int verbose);
#  endif

/* inline functions */
static inline void Spectab_record(Spectab * tab, unsigned alleleCount,
                                  unsigned wgt);
static inline long unsigned Spectab_nObs(const Spectab * tab);

/*
 * Record a single SNP.  Wgt is the number of copies of this
 * SNP. May be 0 or >1 in bootstrap samples.
 */
static inline void Spectab_record(Spectab * tab, unsigned alleleCount,
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

#endif
