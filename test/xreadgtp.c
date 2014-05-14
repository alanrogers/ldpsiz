/**
 * @file xreadgtp.c
 * @author Alan R. Rogers
 * @brief Test readgtp.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "readgtp.h"
#include "assign.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstInput = "# source              = test input\n\
# Haploid sample size   = 10\n\
# Ploidy                = %d\n\
#   snp_id     nucpos    mappos alleles genotypes\n\
         0        262    0.0262       01 0001\n\
         1        362    0.0362       01 1100\n\
         2        536    0.0536       01 1000\n\
         3        799    0.0799       01 0010\n\
         4        861    0.0861       01 0010\n\
         5       1337    0.1337      01 0110\n\
         6       1564    0.1564      01 1110\n\
         7       1905    0.1905      01 0010\n\
         8       1968    0.1968      01 1001\n\
         9       2419    0.2419      01 0010\n";

#undef NDEBUG
int main(void) {
    long        j;
    int         ploidy = 1;
    const char *tstFname = "readgtp-tmp.gtp";

    FILE       *fp = fopen(tstFname, "w");

    fputs(tstInput, fp);
    fclose(fp);

    FILE       *ifp = fopen(tstFname, "r");

    if(ifp == NULL) {
        fprintf(stderr, "Can't open file %s@%s:%d\n",
                tstFname, __FILE__, __LINE__);
        exit(1);
    }

    /******** test Gtp_readSNP *********/

    double      mappos;
    int         chromosome = 12;
    int         twoNsmp = 99;
    unsigned char sitedat[500];
    char        simcmd[200];

    rewind(ifp);

    Assignment *a = Gtp_readHdr(ifp);

    if(a == NULL) {
        fprintf(stderr, "Gtp_readHdr: no assignments found in header\n");
        unitTstResult("Gtp_readHdr", "FAIL");
    }
    Assignment_setInt(a, "chromosome", &chromosome, !MANDATORY);
    Assignment_setInt(a, "Haploid sample size", &twoNsmp, MANDATORY);
    Assignment_setString(a, "macs cmd", simcmd, sizeof(simcmd), !MANDATORY);
    Assignment_free(a);
    a = NULL;

    long        nsnps = 0;

    j = 0;
    ploidy = 1;                 /* haploid */
    char        snpId[20], alleles[20];

    while(1) {
        j = Gtp_readSNP(ifp, snpId, sizeof(snpId),
                        &mappos,
                        alleles, sizeof(alleles),
                        sitedat, sizeof(sitedat), ploidy == 2);
        if(j == EOF)
            break;

        assert(j = 4);
        assert(nsnps == strtol(snpId, NULL, 10));
        assert(strncmp(alleles, "01", 2) == 0);
        if(nsnps == 9) {
            assert(sitedat[0] == 0);
            assert(sitedat[1] == 0);
            assert(sitedat[2] == 1);
            assert(sitedat[3] == 0);
        }
        ++nsnps;
    }
    assert(nsnps == 10);

    rewind(ifp);

    /* test for Diploid data */
    j = nsnps = 0;
    ploidy = 2;
    while(1) {
        j = Gtp_readSNP(ifp, snpId, sizeof(snpId),
                        &mappos,
                        alleles, sizeof(alleles),
                        sitedat, sizeof(sitedat), ploidy == 2);
        if(j == EOF)
            break;

        assert(j = 2);
        assert(nsnps == strtol(snpId, NULL, 10));
        assert(strncmp(alleles, "01", 2) == 0);
        if(nsnps == 9) {
            assert(sitedat[0] == 0);
            assert(sitedat[1] == 2);
        }
        ++nsnps;
    }
    assert(nsnps == 10);
    fclose(ifp);
    unitTstResult("readgtp", "OK");

    remove(tstFname);           /* delete file from directory */

    return 0;
}
