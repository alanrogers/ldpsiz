/**
 * @file readgtp.c
 * @author Alan R. Rogers
 * @brief Parse files written in .gtp format, which are used as input by obsld.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include <limits.h>
#include "tokenizer.h"
#include "misc.h"
#include "fileindex.h"
#include "assign.h"
#include "readgtp.h"

/**
 * Read header of file ifp in order to parse assignment
 * statements. These assignments are stored an a structure of class
 * Assignment, which is returned to the calling function.
 *
 * @param[in] ifp Points to input file stream.
 * @returns pointer to Assignment structure.
 * @sideeffect The file pointer ifp is positioned so that the next
 * line read will be the first SNP.
 */
Assignment *Gtp_readHdr(FILE * ifp) {
    char        buff[500];
    int         ntokens = -1;
    Tokenizer  *tkz = Tokenizer_new(50);
    Assignment *a = NULL;

    while(1) {
        if(fgets(buff, sizeof(buff), ifp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        if(strstr(buff, "snp_id") != NULL)  /* end of header */
            break;

        Tokenizer_split(tkz, buff, "=");    /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens != 2)
            continue;           /* not an assignment, so skip this line */

        /*
         * Read assignment statements, store in Assignment structure.
         */
        a = Assignment_new(a,
                           Tokenizer_token(tkz, 0), Tokenizer_token(tkz, 1));
    }

    Tokenizer_free(tkz);
    return a;
}

/**
 * Read a SNP from the GTP output file. Return 0 on success, 1 on failure.
 *
 * @param[in] ifp Points to input file stream.
 * @param[out] snpId Character buffer to hold the string that
 * identifies current SNP.
 * @param[in] snpIdSize Size of snpId. If snpIdSize==0, nothing will
 * be written into snpId. In this case, it is OK to pass NULL for the
 * value of snpId.
 * @param[out] mappos Map position of current SNP in centimorgans.
 * @param[out] alleles This string will contain a list of the alleles
 * present at current SNP.
 * @param[in] allelesSize Size of "alleles". If allelesSize==0, nothing will
 * be written into alleles. In this case, it is OK to pass NULL for the
 * value of alleles.
 * @param[out] sitedat Unsigned character buffer whose i'th entry will
 * hold the genotype data for the i'th individual. This value is an integer
 * that encodes genotypic values as described in the documentation of
 * function encodeGtype, defined in file misc.h.
 * @param[in] sitedatSize Size of the sitedat array.
 * @param[in] isDiploid The data are assumed to be diploid if
 * isDiploid is nonzero. Otherwise, haploidy is assumed.
 * @returns A positive number (the number of genotypes) on success,
 * but a negative number (EOF) on end of file.
 */
int Gtp_readSNP(FILE * ifp,
                char *snpId, size_t snpIdSize,
                double *mappos,
                char *alleles, size_t allelesSize,
                unsigned char *sitedat, size_t sitedatSize, int isDiploid) {
    int         ntokens = 0;
    unsigned    nalleles;
    unsigned    ngtypes = 0;
    char        buff[500];

    /*
     * Might be worth allocating this in calling program to avoid
     * repeated calls to malloc.
     */
    Tokenizer  *tkz = Tokenizer_new(50);

    /*
     * Repeat this loop until we find an input line with 4
     * tokens. Assume that such a line is a valid SNP.  Then count the
     * number of alleles in the genotype string. If ntokens != 5, the
     * line cannot be a SNP, so keep reading.  If the number of
     * alleles is less than 2, keep reading.
     */
    do {
        if(fgets(buff, sizeof(buff), ifp) == NULL) {
            Tokenizer_free(tkz);
            return EOF;
        }

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        Tokenizer_split(tkz, buff, " \t");
        ntokens = Tokenizer_strip(tkz, " \t\n");

        nalleles = 0;
        if(ntokens == 5) {
            if(strchr(Tokenizer_token(tkz, 4), '1'))
                ++nalleles;
            if(strchr(Tokenizer_token(tkz, 4), '0'))
                ++nalleles;
            if(strchr(Tokenizer_token(tkz, 4), 'h'))
                nalleles = 2;
        }
    } while(ntokens != 5 || nalleles < 2);

    if(snpId)
        snprintf(snpId, snpIdSize, "%s", Tokenizer_token(tkz, 0));
    if(mappos)
        *mappos = strtod(Tokenizer_token(tkz, 2), NULL);
    if(alleles)
        snprintf(alleles, allelesSize, "%s", Tokenizer_token(tkz, 3));
    if(sitedat) {
        if(isDiploid)
            ngtypes = encodeDiploid(sitedat, sitedatSize,
                                    Tokenizer_token(tkz, 4));
        else
            ngtypes = encodeHaploid(sitedat, sitedatSize,
                                    Tokenizer_token(tkz, 4));
    }

    /* snprintf(sitedat, sitedatSize, "%s", Tokenizer_token(tkz, 3)); */

    Tokenizer_free(tkz);
    return (int) ngtypes;
}
