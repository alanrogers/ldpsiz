/**
 * @file boot.h
 * @author Alan R. Rogers
 * @brief Header for boot.c.
 * @copyright Copyright (c) 2014, 2015 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_BOOT_H
#define LDPSIZ_BOOT_H

#include "typedefs.h"
#include <gsl/gsl_rng.h>
void        confidenceBounds(double *lowBnd, double *highBnd,
                             double confidence, double *v, long len);
double      interpolate(double p, double *v, long len);
Boot       *Boot_new(long nSNPs, long nReps, unsigned twoNsamp,
                     int folded, long blockLength,
                     double windowcm, int nBins, gsl_rng * rng);
void        Boot_addLD(Boot * boot, double Dsq, double pqpq, double sep_cm,
                       const SNP * snp1, const SNP * snp2);
void        Boot_addAlleleCount(Boot * boot, unsigned x, const SNP * snp);
void        Boot_free(Boot * boot);
int         Boot_equals(const Boot * x, const Boot * y);
Boot       *Boot_dup(const Boot * old);
long        Boot_nReps(const Boot * boot);
int         Boot_nBins(const Boot * boot);
long        Boot_nBlocks(const Boot * boot);
long        Boot_nSNPs(const Boot * boot);
void        Boot_plus_equals(Boot * x, const Boot * y);
void        Boot_dump(const Boot * boot, FILE * ofp);
Boot       *Boot_restore(FILE * ifp);
long        Boot_multiplicity(const Boot * boot, long ndx, long rep);
void        Boot_get_rep(Boot * boot, DblArray *sigdsq, DblArray *rsq,
                         DblArray *cm, ULIntArray *nobs,
                         ULIntArray *spectrum, int rep);
long unsigned Boot_rawCounts(const Boot * boot, int rep, int bin,
                             double *numerator, double *denominator,
                             double *sumRsq, double *sep_cm);
long        Boot_purge(Boot * boot);
void        Boot_print(const Boot * boot, FILE * ofp);

#ifndef NDEBUG
unsigned Boot_multiplicity_slow(Boot * boot, long snp, long rep);
#endif

BootConf   *BootConf_new(Boot * boot, double confidence);
void        BootConf_printHdr(const BootConf * bc, FILE * ofp);
double      BootConf_lowBound(const BootConf * bc, long bin);
double      BootConf_highBound(const BootConf * bc, long bin);
double      BootConf_loSpecBound(const BootConf * bc, long i);
double      BootConf_hiSpecBound(const BootConf * bc, long k);
void        BootConf_print(const BootConf * bc, FILE * ofp);
void        BootConf_free(BootConf * bc);

#endif
