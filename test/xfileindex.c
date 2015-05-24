/**
 * @file xfileindex.c
 * @author Alan R. Rogers
 * @brief Test fileindex.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "fileindex.h"
#include "misc.h"
#include "readgtp.h"
#include "assign.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

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

int main(int argc, char **argv) {
    int         rval, verbose = 0;
    int         bad, ploidy = 1;
    long        i, j;
    FileIndex  *fndx = FileIndex_new();
    SNPLoc     *slp;
    char        snpId[20];
    unsigned char sitedat[1000];
    const unsigned nptrs = FileIndex_nPtrs();
    const unsigned chunksize = FileIndex_chunkSize();

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xfileindex [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xfileindex [-v]\n");
    }

    const char *tstFname = "fileindex-tmp.gtp";

    FILE       *fp = fopen(tstFname, "w");

    fprintf(fp, tstInput, ploidy);
    fclose(fp);

    FILE       *ifp = fopen(tstFname, "r");

    if(ifp == NULL) {
        fprintf(stderr, "Can't open file %s@%s:%d\n",
                tstFname, __FILE__, __LINE__);
        exit(1);
    }

    assert(FileIndex_nSNPs(fndx) == 0);

    FileIndex_realloc(fndx, nptrs - 1);
    assert(FileIndex_nPtrsAllocated(fndx) == nptrs);
    i = 100 + nptrs;
    FileIndex_realloc(fndx, i);
    j = (i / nptrs) * nptrs;
    if(j < i)
        j += nptrs;
    assert(FileIndex_nPtrsAllocated(fndx) == j);

    FileIndex_free(fndx);
    fndx = FileIndex_new();

	/* this loop takes a lot of time */
    for(i = 0; i < (nptrs * chunksize + 1); ++i)
        FileIndex_push(fndx, i, i);

	/* sanity check is slow too */
    FileIndex_sanityCheck(fndx, __FILE__, __LINE__);

    bad = 0;
    for(i = 0; i < (nptrs * chunksize + 1); ++i) {
        slp = FileIndex_atNdx(fndx, i);
		long seekpos = SNPLoc_seekpos(slp);
		double mappos = SNPLoc_mappos(slp);
        if(i != seekpos) {
            printf("ERR: i=%ld seekpos=%ld\n", i, seekpos);
            bad = 1;
        }
        if(i != mappos) {
            printf("ERR: i=%ld mappos=%lg\n", i, mappos);
            bad = 1;
        }
    }
    if(bad) {
        unitTstResult("FileIndex", "FAIL");
        die("seekpos and index values don't match", __FILE__, __LINE__);
    }

    bad = 0;
    for(i = 0; i < chunksize + 1; ++i) {
        double      mappos = FileIndex_getMapPos(fndx, i);

        if(i != FileIndex_atMapPos(fndx, mappos)) {
            fprintf(stderr, "%s:%d: FileIndex_atMapPos doesn't invert"
                    "FileIndex_getMapPos\n", __FILE__, __LINE__);
            bad = 1;
        }

        j = FileIndex_getSeekPos(fndx, i);
        if(j != SNPLoc_seekpos(FileIndex_atNdx(fndx, i))) {
            fprintf(stderr, "%s:%d: FileIndex_atSeekPos doesn't match"
                    "FileIndex_atNdx\n", __FILE__, __LINE__);
            bad = 1;
        }

    }
    if(bad) {
        unitTstResult("FileIndex", "FAIL");
        die("aborting because of errors", __FILE__, __LINE__);
    }

    /******** test FileIndex_readFile *********/

    rewind(ifp);
    FileIndex_free(fndx);
    fndx = FileIndex_readFile(ifp);

    if(fndx == NULL) {
        printf("  NULL return from FileIndex_readFile\n");
        unitTstResult("FileIndex_readFile", "FAIL");
    } else {
        long        ndx1;
        double      mappos1, mappos2;
        SNPLoc     *sl;

        FileIndex_sanityCheck(fndx, __FILE__, __LINE__);

        if(FileIndex_nSNPs(fndx) != 10) {
            printf("FileIndex_nSNPs = %ld; should be 10\n",
                   FileIndex_nSNPs(fndx));
            FileIndex_print(fndx, stdout);
            unitTstResult("FileIndex_readFile", "FAIL");
            fflush(stdout);
        }
        assert(FileIndex_nSNPs(fndx) == 10);

        /* check 0th SNP */
        ndx1 = 0;
        sl = FileIndex_atNdx(fndx, ndx1);
        mappos1 = SNPLoc_mappos(sl);   // as recorded in fndx
        fseek(ifp, SNPLoc_seekpos(sl), SEEK_SET);
        rval = Gtp_readSNP(ifp, snpId, sizeof(snpId), &mappos2, NULL,
                           0,    // ignore list of alleles
                           sitedat, sizeof(sitedat), ploidy == 2);
        assert(rval > 0);
        assert(ndx1 == strtol(snpId, NULL, 10));
        assert(mappos1 == mappos2);

        /* check last SNP */
        ndx1 = FileIndex_nSNPs(fndx) - 1;
        sl = FileIndex_atNdx(fndx, ndx1);
        mappos1 = SNPLoc_mappos(sl);   //* as recorded in fndx
        fseek(ifp, SNPLoc_seekpos(sl), SEEK_SET);
        rval = Gtp_readSNP(ifp, snpId, sizeof(snpId), &mappos2, NULL,
                           0,    // ignore list of alleles
                           sitedat, sizeof(sitedat), ploidy == 2);
        assert(rval > 0);
        assert(ndx1 == strtol(snpId, NULL, 10));
        assert(mappos1 == mappos2);

        /* check middle SNP */
        ndx1 /= 2;
        sl = FileIndex_atNdx(fndx, ndx1);
        mappos1 = SNPLoc_mappos(sl);   /* as recorded in fndx */
        fseek(ifp, SNPLoc_seekpos(sl), SEEK_SET);
        rval = Gtp_readSNP(ifp, snpId, sizeof(snpId), &mappos2, NULL,
                           0,    /* ignore list of alleles */
                           sitedat, sizeof(sitedat), ploidy == 2);
        assert(rval > 0);
        assert(ndx1 == strtol(snpId, NULL, 10));
        assert(mappos1 == mappos2);
    }

    assert(ploidy == FileIndex_ploidy(fndx));
    assert(4 == FileIndex_nGtype(fndx));

    if(verbose)
        FileIndex_print(fndx, stdout);
    unitTstResult("FileIndex", "OK");

    int         nthreads = 2;
    double      windowsize = 0.2;

    ThreadBounds *tb = ThreadBounds_new(nthreads, windowsize, fndx);

    if(verbose)
        ThreadBounds_print(tb, nthreads, stdout);
    ThreadBounds_sanityCheck(tb, __FILE__, __LINE__);

    long        ndx[3], seekpos[3];
    double      mappos[3];

    for(i = 0; i < nthreads; ++i) {
        for(j = 0; j < 3; ++j) {
            seekpos[j] = tb[i].seekpos[j];
            ndx[j] = tb[i].ndx[j];
            fseek(ifp, seekpos[j], SEEK_SET);
            rval = Gtp_readSNP(ifp, snpId, sizeof(snpId), mappos + j, NULL,
                               0,  // ignore list of alleles
                               sitedat, sizeof(sitedat), ploidy == 2);
            assert(mappos[j] == tb[i].mappos[j]);
            if(j > 0) {
                assert(ndx[j] >= ndx[j - 1]);
                assert(seekpos[j] >= seekpos[j - 1]);
                assert(mappos[j] >= mappos[j - 1]);
            } else {
                assert(ndx[j] >= 0);
                assert(seekpos[j] >= 0);
                assert(mappos[j] >= 0.0);
            }
            switch (j) {
            case 0:
                assert(ndx[j] == ThreadBounds_ndx_initial(tb + i));
                assert(seekpos[j] == ThreadBounds_seekpos_initial(tb + i));
                break;
            case 1:
                assert(ndx[j] == ThreadBounds_ndx_1stFocal(tb + i));
                assert(seekpos[j] == ThreadBounds_seekpos_1stFocal(tb + i));
                break;
            case 2:
                assert(ndx[j] == ThreadBounds_ndx_lastFocal(tb + i));
                assert(seekpos[j] == ThreadBounds_seekpos_lastFocal(tb + i));
                break;
            }
        }

    }

    unitTstResult("ThreadBounds", "OK");

    FileIndex_free(fndx);
    ThreadBounds_free(tb);
    fclose(ifp);
    remove(tstFname);           /* delete file from directory */

    return 0;
}
