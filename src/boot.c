/**
 * @file boot.c
 * @author Alan R. Rogers
 * @brief Functions for a moving blocks bootstrap.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "boot.h"
#include "tabulation.h"
#include "spectab.h"
#include "misc.h"
#include "window.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

/**
 * Contains the all data involved in a moving blocks bootstrap.
 */
struct Boot {
    long        nSNPs;       // number of polymorphic sites in data
    long        nReps;       // number of bootstrap replicates 
    long        blockLength; // number of SNPs per block
    long        nBlocks;     // number of blocks in each replicate
    int         nBins;       // tabulate recombination rates into this
                             // number of bins 
    long      **start;       // start[i][j] = start of j'th block in i'th rep
    Tabulation **tab;        // LD tabulation for each replicate
    Spectab    **spectab;    // spectrum tabulation for each replicate
};

/** Contains the data for a bootstrap confidence interval. */
struct BootConf {
    long        nReps;       // repetitions
    long        blockLength; // nucleotide positions per block
    int         nBins;       // # of recombination rate bins
    unsigned    specDim;     // dimension of spectrum
    double      confidence;  // size of confidence region
    double     *low, *high;  // LD confidence bounds
    double     *loSpec, *hiSpec; // spectrum confidence bounds
    double     *sep_cm;      // sep_cm[i]=mean centimorgans separating
                             // pairs of sites within bin i. 
};

/**
 * Allocate Boot's arrays. This code is used in several places, and I
 * have it here to ensure consistency.
 */
static void Boot_allocArrays(Boot * boot) {

    long        i;

    boot->start = calloc((unsigned long) boot->nReps, sizeof(boot->start[0]));
    checkmem(boot->start, __FILE__, __LINE__);

    boot->tab = calloc((unsigned long) boot->nReps, sizeof(boot->tab[0]));
    checkmem(boot->tab, __FILE__, __LINE__);

    boot->spectab = calloc((unsigned long) boot->nReps,
                           sizeof(boot->spectab[0]));
    checkmem(boot->spectab, __FILE__, __LINE__);

    for(i = 0; i < boot->nReps; ++i) {
        boot->start[i] = calloc((unsigned long) boot->nBlocks,
                                sizeof(boot->start[0][0]));
        checkmem(boot->start[i], __FILE__, __LINE__);
    }
}

/// Constructor for class Boot.
Boot       *Boot_new(long nSNPs, long nReps, unsigned twoNsamp,
                     int folded, long blockLength,
                     double windowcm, int nBins, gsl_rng * rng) {

    if(nReps==0)
        return NULL;

    if(blockLength > nSNPs) {
        fprintf(stderr,
                "ERR@%s:%d: in Boot_new, nSNPs must be >blockLength.\n"
                " Instead, nSNPs=%ld, blockLength=%ld.\n",
                __FILE__, __LINE__, nSNPs, blockLength);
        fprintf(stderr, " Use --blocksize argument to reduce blockLength.\n");
        exit(1);
    }

    long        i, j;

    /* Block positions are uniform on [0, nSNPs-blockLength+1). */
    unsigned long endpos = nSNPs - blockLength + 1;

    Boot       *boot = malloc(sizeof(Boot));
    checkmem(boot, __FILE__, __LINE__);

    boot->nSNPs = nSNPs;
    boot->nReps = nReps;
    boot->blockLength = blockLength;
    boot->nBlocks = (long) round(nSNPs / ((double) blockLength));
    assert(boot->nBlocks > 0);
    boot->nBins = nBins;

    Boot_allocArrays(boot);

    for(i = 0; i < boot->nReps; ++i) {
        boot->tab[i] = Tabulation_new(windowcm, nBins);
        checkmem(boot->tab[i], __FILE__, __LINE__);

        boot->spectab[i] = Spectab_new(twoNsamp, folded);
        checkmem(boot->spectab[i], __FILE__, __LINE__);

        for(j = 0; j < boot->nBlocks; ++j)
            boot->start[i][j] = gsl_rng_uniform_int(rng, endpos);

        qsort(boot->start[i], (size_t) boot->nBlocks,
              sizeof(boot->start[0][0]), compareLongs);
    }
    return boot;
}

/*
 * How many copies of snp are present in a given repetition (rep)?
 */
long Boot_multiplicity(const Boot * boot, long snp, long rep) {
    long        lndx, hndx, lowtarget;

    myassert(snp < boot->nSNPs);

    /* lndx is index of first block containing snp */
    lowtarget = snp - boot->blockLength + 1;
    lndx = long_first_geq(lowtarget, boot->start[rep], boot->nBlocks);
    if(lndx == boot->nBlocks || boot->start[rep][lndx] > snp)
        return 0;

    myassert(snp >= boot->start[rep][lndx]);
    myassert(snp - boot->start[rep][lndx] < boot->blockLength);

    /* hndx is index of first block not containing snp */
    hndx = long_first_geq(snp + 1, boot->start[rep] + lndx,
                          boot->nBlocks - lndx);
    hndx += lndx;

    myassert(hndx == 0
             || boot->start[rep][hndx - 1] - snp < boot->blockLength);

    return hndx - lndx;
}

/*
 * Add one LD value to a Boot structure. On entry, boot points to the
 * Boot structure, Dsq and pqpq are the contributions to the numerator
 * and denominator of sigma_d^2, sep_cm is the separation in kilobases
 * between the two loci, and ndx1 and ndx2 are their positions within
 * the list of SNPs.
 */
void Boot_addLD(Boot * boot, double Dsq, double pqpq, double sep_cm,
                const SNP * snp1, const SNP * snp2) {
    for(register int rep = 0; rep < boot->nReps; ++rep) {
        register unsigned wgt = SNP_multiplicity(snp1, rep)
            * SNP_multiplicity(snp2, rep);

        if(wgt == 0)
            continue;
        Tabulation_record(boot->tab[rep], Dsq, pqpq, sep_cm, wgt);
    }
}

/**
 * Add one allele count to a Boot structure. On entry, boot points to the
 * Boot structure, x is the number of copies of an allele at the SNP
 * whose position is given by ndx.
 */
void Boot_addAlleleCount(Boot * boot, unsigned x, const SNP * snp) {
    for(register int rep = 0; rep < boot->nReps; ++rep) {
        register unsigned wgt = SNP_multiplicity(snp, rep);
        if(wgt == 0)
            continue;
        Spectab_record(boot->spectab[rep], x, wgt);
    }
}

/*
 * Allocate and initialize a duplicate of a Boot structure. Deallocate
 * with Boot_free.
 */
Boot       *Boot_dup(const Boot * old) {
    myassert(old != NULL);
    myassert(old->nBlocks > 0);
    Boot       *new = memdup(old, sizeof(*old));

    myassert(new->nBlocks > 0);

    new->start = malloc(new->nReps * sizeof(new->start[0]));
    checkmem(new->start, __FILE__, __LINE__);

    new->tab = malloc(new->nReps * sizeof(new->tab[0]));
    checkmem(new->tab, __FILE__, __LINE__);

    new->spectab = malloc(new->nReps * sizeof(new->spectab[0]));
    checkmem(new->spectab, __FILE__, __LINE__);

    myassert(sizeof(new->start[0][0]) > 0);

    for(int i = 0; i < new->nReps; ++i) {
        new->tab[i] = Tabulation_dup(old->tab[i]);
        new->spectab[i] = Spectab_dup(old->spectab[i]);
        new->start[i] = memdup(old->start[i],
                               new->nBlocks * sizeof(new->start[0][0]));
    }
    return new;
}

/* Return number of bootstrap repetitions */
long Boot_nReps(const Boot * boot) {
    myassert(boot);

    return boot->nReps;
}

/* Return number of bins */
int Boot_nBins(const Boot * boot) {
    myassert(boot);

    return boot->nBins;
}

/** Return number of blocks */
long Boot_nBlocks(const Boot * boot) {
    myassert(boot);

    return boot->nBlocks;
}

/** Return number of SNPs */
long Boot_nSNPs(const Boot * boot) {
    myassert(boot);

    return boot->nSNPs;
}

/* deallocate bootstrap */
void Boot_free(Boot * boot) {
    long        i;

    myassert(boot != NULL);

    for(i = 0; i < boot->nReps; ++i) {
        free(boot->start[i]);
        Tabulation_free(boot->tab[i]);
        Spectab_free(boot->spectab[i]);
        boot->start[i] = NULL;
        boot->tab[i] = NULL;
    }
    free(boot->start);
    free(boot->tab);
    free(boot->spectab);
    boot->start = NULL;
    boot->tab = NULL;
    free(boot);
}

/**
 * Interpolate in order to approximate the value v[p*(len-1)].
 * Return NaN if len==0.
 */
double interpolate(double p, double *v, long len) {
    if(len == 0)
        return strtod("NAN", 0);
    long        i, j;
    double      w;
    double      goal = p * (len - 1);

    i = floor(goal);
    j = ceil(goal);

    myassert(i >= 0);
    myassert(j < len);

    if(i == j)                  /* no interpolation needed */
        return v[i];
    w = goal - i;
    return (1.0 - w) * v[i] + w * v[j];
}

/* Return 1 if the two Boot structs are idential; zero otherwise */
int Boot_equals(const Boot * x, const Boot * y) {
    if(x == NULL || y == NULL)
        return false;
    if(x->nSNPs != y->nSNPs)
        return false;
    if(x->nReps != y->nReps)
        return false;
    if(x->blockLength != y->blockLength)
        return false;
    if(x->nBins != y->nBins)
        return false;

    long        i, j;

    for(i = 0; i < x->nReps; ++i) {
        if(!Tabulation_equals(x->tab[i], y->tab[i]))
            return false;

        if(!Spectab_equals(x->spectab[i], y->spectab[i]))
            return false;

        for(j = 0; j < x->nBlocks; ++j)
            if(x->start[i][j] != y->start[i][j])
                return false;
    }
    return true;
}

/*
 * Add the contents of Boot structure y to Boot structure x.
 * On return x will summarize all the LD values that were originally
 * contained in the two original Boot structures.
 */
void Boot_plus_equals(Boot * x, const Boot * y) {
    int         rep;

    if(x == NULL)
        die("Boot_plus_equals: destination is NULL", __FILE__, __LINE__);
    if(y == NULL)
        return;

    if(x->nSNPs != y->nSNPs
       || x->nReps != y->nReps
       || x->blockLength != y->blockLength
       || x->nBlocks != y->nBlocks || x->nBins != y->nBins)
        die("Boot_plus_equals: unconformable arguments", __FILE__, __LINE__);

    for(rep = 0; rep < x->nReps; ++rep) {
        long        block;

        for(block = 0; block < x->nBlocks; ++block) {
            if(x->start[rep][block] != y->start[rep][block])
                die("Boot_plus_equals: incompatible block starts",
                    __FILE__, __LINE__);
        }
        Tabulation_plus_equals(x->tab[rep], y->tab[rep]);
        Spectab_plus_equals(x->spectab[rep], y->spectab[rep]);
    }
    return;
}

/*
 * Write a text file representing a Boot structure. To restore the
 * Boot structure from this file, use Boot_restore.
 */
void Boot_dump(const Boot * boot, FILE * ofp) {
    long        rep, block;

    if(fprintf(ofp, " %ld %ld %ld %ld %d\n",
               boot->nSNPs, boot->nReps, boot->blockLength,
               boot->nBlocks, boot->nBins) < 0)
        eprintf("fprintf", __FILE__, __LINE__);

    for(rep = 0; rep < boot->nReps; ++rep) {
        for(block = 0; block < boot->nBlocks; ++block) {
            if(fprintf(ofp, " %ld", boot->start[rep][block]) < 0)
                eprintf("fprintf", __FILE__, __LINE__);
        }
        putc('\n', ofp);
        Tabulation_dump(boot->tab[rep], ofp);
        Spectab_dump(boot->spectab[rep], ofp);
    }
}

/*
 * Create a Boot structure by reading a file produced by
 * Boot_dump.
 */
Boot       *Boot_restore(FILE * ifp) {

    long        nSNPs, nReps, blockLength, nBlocks, rep, block;
    int         nBins, rval;
    Boot       *boot = malloc(sizeof(Boot));

    checkmem(boot, __FILE__, __LINE__);

    rval = fscanf(ifp, "%ld %ld %ld %ld %d",
                  &nSNPs, &nReps, &blockLength, &nBlocks, &nBins);
    if(rval != 5)
        eprintf("ERR@%s:%d: fscanf returned %d instead of 5\n",
                __FILE__, __LINE__, rval);

    boot->nSNPs = nSNPs;
    boot->nReps = nReps;
    boot->blockLength = blockLength;
    boot->nBlocks = (long) round(nSNPs / ((double) blockLength));
    boot->nBins = nBins;

    if(nBlocks != boot->nBlocks) {
        fprintf(stderr, "from file:"
                " nSNPs=%ld nReps=%ld blockLength=%ld nBlocks=%ld nBins=%d\n",
                nSNPs, nReps, blockLength, nBlocks, nBins);
        fprintf(stderr, "in boot  :"
                " nSNPs=%ld nReps=%ld blockLength=%ld nBlocks=%ld nBins=%d\n",
                boot->nSNPs, boot->nReps, boot->blockLength,
                boot->nBlocks, boot->nBins);
        eprintf("ERR@%s:%d: nBlocks=%ld but boot->nBlocks=%ld\n",
                __FILE__, __LINE__, nBlocks, boot->nBlocks);
    }

    Boot_allocArrays(boot);

    for(rep = 0; rep < boot->nReps; ++rep) {
        for(block = 0; block < boot->nBlocks; ++block) {
            rval = fscanf(ifp, " %ld", boot->start[rep] + block);
            if(rval != 1)
                eprintf("ERR@%s:%d: fscanf returned %d instead of 1\n",
                        __FILE__, __LINE__, rval);
        }
        boot->tab[rep] = Tabulation_restore(ifp);
        boot->spectab[rep] = Spectab_restore(ifp);
        if(!Tabulation_isfinite(boot->tab[rep])) {
            fprintf(stderr,
                    "Warning@%s:%d: bootstrap replicate %ld is non-finite\n",
                    __FILE__, __LINE__, rep);
        }
    }
    return boot;
}

/**
 * Fill arrays sigdsq, cm, and nobs with values for bootstrap
 * repetition "rep". If nobs==NULL, nothing is stored there.
 */
void Boot_get_rep(Boot * boot, DblArray *sigdsq, DblArray *rsq,
                  DblArray *cm, ULIntArray *nobs,
                  ULIntArray *spectrum, int rep) {
    myassert(boot);
    myassert(sigdsq);
    myassert(cm);
    myassert(rep < boot->nReps);
    myassert(rep >= 0);
    Tabulation_report(boot->tab[rep], cm, nobs, sigdsq, rsq);
    Spectab_report(boot->spectab[rep], spectrum);
}

/*
 * Return the raw counts for a single rep and bin. Counts for
 * numerator, denominator, and sep_cm are returned in the
 * corresponding arguments. The function itself returns nobs.
 */
long unsigned Boot_rawCounts(const Boot * boot, int rep, int bin,
                             double *numerator, double *denominator,
                             double *sumRsq, double *sep_cm) {

    assert(rep < boot->nReps);
    assert(bin < boot->nBins);

    return Tabulation_rawCounts(boot->tab[rep], bin, numerator,
                                denominator, sumRsq, sep_cm);
}

/**
 * Remove replicates in which some bins have zero observations.
 * These would generate NaN values in the calculation of sigdsq,
 * which are not handled by the minimizer. Return revised number of
 * bootstrap replicates.
 */
long Boot_purge(Boot * boot) {
    myassert(boot);

    long        rep, nGoodReps = boot->nReps;
    int         bin;

    rep = 0;
    while(rep < nGoodReps) {
        int         clean = 1;

        for(bin = 0; bin < boot->nBins; ++bin) {
            if(Tabulation_nObs(boot->tab[rep], bin) == 0) {
                /*
                 * Current rep has an empty bin, which is not OK.
                 * Mark this rep as dirty and break out of loop.
                 */
                clean = 0;
                break;
            }
            if(Tabulation_overflow(boot->tab[rep])) {
                /*
                 * Current rep has invalid data because its Tabulation
                 * overflowed.
                 */
                clean = 0;
                break;
            }
        }
        /*
         * If current rep is clean, then increment reps; otherwise
         * free current rep, move in terminal rep, and reduce count of
         * good reps.
         */
        if(clean)
            ++rep;
        else {
            free(boot->start[rep]);
            Tabulation_free(boot->tab[rep]);
            Spectab_free(boot->spectab[rep]);

            if(rep < nGoodReps - 1) {
                boot->start[rep] = boot->start[nGoodReps - 1];
                boot->tab[rep] = boot->tab[nGoodReps - 1];
            }
            --nGoodReps;
        }
    }
    boot->nReps = nGoodReps;
    return boot->nReps;
}

/** Print a Boot object */
void Boot_print(const Boot * boot, FILE * ofp) {
    long        rep, block;

    fprintf(ofp,
            "Boot_print: nSNPs=%ld nReps=%ld blockLength=%ld nBlocks=%ld\n",
            boot->nSNPs, boot->nReps, boot->blockLength, boot->nBlocks);

    fprintf(ofp, "Block starts:\n");
    for(rep = 0; rep < boot->nReps; ++rep) {
        fprintf(ofp, "  rep %ld:", rep);
        for(block = 0; block < boot->nBlocks; ++block)
            fprintf(ofp, " %ld", boot->start[rep][block]);
        putc('\n', ofp);
        // Tabulation_print(boot->tab[rep], ofp);
        // Spectab_print(boot->spectab[rep], ofp);
    }
}

#ifndef NDEBUG
/* For debugging Boot_multiplicity */
unsigned Boot_multiplicity_slow(Boot * boot, long snp, long rep) {
    unsigned        i, n = 0;

    for(i = 0; i < boot->nBlocks; ++i) {
        long        distance = snp - boot->start[rep][i];

        if(distance < 0)
            break;
        if(distance < boot->blockLength)
            ++n;
    }
    return n;
}
#endif

BootConf   *BootConf_new(Boot * boot, double confidence) {
    BootConf   *bc = malloc(sizeof(BootConf));

    checkmem(bc, __FILE__, __LINE__);

    bc->confidence = confidence;
    bc->nReps = boot->nReps;
    bc->blockLength = boot->blockLength;
    bc->nBins = boot->nBins;

    bc->low = malloc(bc->nBins * sizeof(bc->low[0]));
    checkmem(bc->low, __FILE__, __LINE__);

    bc->high = malloc(bc->nBins * sizeof(bc->high[0]));
    checkmem(bc->high, __FILE__, __LINE__);

    bc->sep_cm = malloc(bc->nBins * sizeof(bc->sep_cm[0]));
    checkmem(bc->sep_cm, __FILE__, __LINE__);
    memset(bc->sep_cm, 0, bc->nBins * sizeof(bc->sep_cm[0]));

    assert(boot->spectab && boot->spectab[0]);

    // Dimension of site frequency spectrum
    bc->specDim = Spectab_dim(boot->spectab[0]);

    bc->loSpec = malloc(bc->specDim * sizeof(bc->loSpec[0]));
    checkmem(bc->loSpec, __FILE__, __LINE__);

    bc->hiSpec = malloc(bc->specDim * sizeof(bc->hiSpec[0]));
    checkmem(bc->hiSpec, __FILE__, __LINE__);

    int         i, rep, nvals;
    double     v[bc->nReps];

    // Confidence bounds on sigdsq
    for(i = 0; i < bc->nBins; ++i) {
        double      tmp1, nobs = 0, sigdsq;
        long unsigned tmp2;
        nvals = 0;

        for(rep = 0; rep < bc->nReps; ++rep) {
            sigdsq = Tabulation_sigdsq(boot->tab[rep], i, &tmp1, &tmp2);
            if(isfinite(sigdsq)) {
                myassert(nvals < bc->nReps);
                v[nvals] = sigdsq;
                bc->sep_cm[i] += tmp1;
                nobs += tmp2;
                ++nvals;
            }
        }
        nobs /= nvals;
        bc->sep_cm[i] /= nvals;
        if(nvals < 10) {
            bc->low[i] = bc->high[i] = strtod("NAN", NULL);
        } else
            confidenceBounds(bc->low + i, bc->high + i,
                             confidence, v, nvals);
    }

    // Confidence bounds on spectrum
    for(i = 0; i < bc->specDim; ++i) {
        long unsigned count;
        nvals = 0;

        for(rep = 0; rep < bc->nReps; ++rep) {
            count = Spectab_get(boot->spectab[rep], i);
            v[nvals++] = (double) count;
        }
        if(nvals < 10) {
            bc->low[i] = bc->high[i] = strtod("NAN", NULL);
        } else
            confidenceBounds(bc->loSpec + i, bc->hiSpec + i,
                             confidence, v, nvals);
    }

    return bc;
}

/**
 * Calculate confidence bounds from a vector of values representing
 * samples drawn from the sampling distribution of some estimator.
 *
 * To calculate the lower bound (*lowBnd), the function calculates the
 * total probability mass in the tails (1 - confidence) and divides
 * this into two equal parts to find p, the probability mass in each
 * tail. It then estimates a value L such that a fraction p of the data
 * values are less than or equal to L. To find this value, the
 * function uses linear interpolation between the sorted list of data
 * values.
 *
 * The upper bound (*highBnd) is calculated in an analogous fashion.
 *
 * @param[out] lowBnd,highBnd Calculated results will be written into
 * these memory locations.
 * @param[in] confidence Fraction of sampling distribution that lies
 * inside the confidence bounds.
 * @param[in] v The vector of values.
 * @param[in] len The number of values inf v.
 * @sideeffect The function sorts the vector v.
 */
void confidenceBounds(double *lowBnd, double *highBnd, double confidence,
                      double *v, long len) {
    double      tailProb = (1.0 - confidence) / 2.0;

    qsort(v, (size_t) len, sizeof(v[0]), compareDoubles);
    *lowBnd = interpolate(tailProb, v, len);
    *highBnd = interpolate(1.0 - tailProb, v, len);
}

void BootConf_printHdr(const BootConf * bc, FILE * ofp) {
    fprintf(ofp, "#%12s: %lg%% confidence bounds"
            " based on moving blocks bootstrap\n",
            "loLD, hiLD", 100.0 * bc->confidence);
    fprintf(ofp, "#%12s: nReps=%ld blockLength=%ld nBins=%d\n",
            "Bootstrap parameters", bc->nReps, bc->blockLength, bc->nBins);
}

double BootConf_lowBound(const BootConf * bc, long bin) {
    return bc->low[bin];
}

double BootConf_highBound(const BootConf * bc, long bin) {
    return bc->high[bin];
}

double BootConf_loSpecBound(const BootConf * bc, long i) {
    return bc->loSpec[i];
}

double BootConf_hiSpecBound(const BootConf * bc, long i) {
    return bc->hiSpec[i];
}

void BootConf_print(const BootConf * bc, FILE * ofp) {
    int         i;

    BootConf_printHdr(bc, ofp);
    fprintf(ofp, "%5s %10s %10s\n", "bin", "loLD", "hiLD");
    for(i = 0; i < bc->nBins; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i, bc->low[i], bc->high[i]);

    putc('\n', ofp);
    fprintf(ofp, "%5s %10s %10s\n", "i", "loSpec", "hiSpec");
    for(i = 0; i < bc->specDim; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i+1, bc->loSpec[i], bc->hiSpec[i]);
}

void BootConf_free(BootConf * bc) {
    free(bc->low);
    free(bc->high);
    free(bc->sep_cm);
    free(bc);
}

