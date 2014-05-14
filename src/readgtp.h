/**
 * @file readgtp.h
 * @author Alan R. Rogers
 * @brief Header for readgtp.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_READGTP
#define LDPSIZ_READGTP

#include <stdio.h>
#include "typedefs.h"

Assignment *Gtp_readHdr(FILE * ifp);
int         Gtp_readSNP(FILE * ifp,
                        char *snpId, size_t snpIdSize,
                        double *mappos,
                        char *alleles, size_t allelesSize,
                        unsigned char *sitedat, size_t sitedatSize,
                        int isDiploid);

static inline int recode_0_1(int c);

/* defined here so compiler can inline */

/* Recode character values '0' and '1' as integers 0 and 1 */
static inline int recode_0_1(int c) {
    return c - '0';
}

#endif                       /* LDPSIZ_READGTP */
