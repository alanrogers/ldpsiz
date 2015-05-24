/**
 * @file fileindex.c
 * @author Alan R. Rogers
 * @brief Class FileIndex indexes of the SNPs in a gtp data file.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "misc.h"
#include "readgtp.h"
#include "assign.h"
#include "fileindex.h"

/** The location of a SNP in a file and on a chromosome */
struct SNPLoc {
    long        seekpos;        /* file position, as reported by ftell */
    double      mappos;         /* map position in centimorgans */
};

/**
 * Index of SNPs in a file. Allows translation between the SNP's
 * index, mappos, and seekpos.
 */
struct FileIndex {
    /*
     * snploc       points to an array of nptrs*sizeof(SNPLoc *)
     * snploc[i]    may be NULL (if not allocated) or points to an array
     *              of CHUNKSIZE*sizeof(SNPLoc)
     * snploc[i][j] is an object of type SNPLoc, not a pointer.
     */
    long        chunksAllocated, ptrsAllocated, nSNPs;
    SNPLoc    **snploc;         /* snploc[i][j] is j'th snp w/i i'th chunk */
    unsigned    ploidy, nGtype;
};

/* number of objects allocated by each malloc call within FileIndex */
#define CHUNKSIZE 4000
#define NPTRS 4000

void        SNPLoc_print(const SNPLoc * sl, FILE * ofp);
void        SNPLoc_sanityCheck(const SNPLoc * sl, const char *file,
                               int lineno);

unsigned FileIndex_nPtrs(void) {
    return NPTRS;
}

unsigned FileIndex_chunkSize(void) {
    return CHUNKSIZE;
}

long FileIndex_nPtrsAllocated(const FileIndex *fndx) {
    return fndx->ptrsAllocated;
}

void SNPLoc_print(const SNPLoc * sl, FILE * ofp) {
    fprintf(ofp, "seekpos=%ld mappos=%lf\n", sl->seekpos, sl->mappos);
}

void SNPLoc_sanityCheck(const SNPLoc * sl, const char *file, int lineno) {
    REQUIRE(sl != NULL, file, lineno);
    REQUIRE(sl->seekpos >= 0, file, lineno);
    REQUIRE(sl->mappos >= 0.0, file, lineno);
}

long SNPLoc_seekpos(const SNPLoc *sl) {
	return sl->seekpos;
}

double SNPLoc_mappos(const SNPLoc *sl) {
	return sl->mappos;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
void FileIndex_sanityCheck(FileIndex * fndx, const char *file, int lineno) {

    REQUIRE(fndx != NULL, file, lineno);

    long        i, nSNPs = FileIndex_nSNPs(fndx);

    REQUIRE(nSNPs >= 0, file, lineno);

    if(nSNPs == 0)
        return;

    SNPLoc     *snp;
    SNPLoc     *snp0 = FileIndex_atNdx(fndx, 0);

    SNPLoc_sanityCheck(snp0, __FILE__, __LINE__);

    long        seekpos = FileIndex_getSeekPos(fndx, 0);

    REQUIRE(seekpos == snp0->seekpos, file, lineno);

    double      mappos = FileIndex_getMapPos(fndx, 0);

    REQUIRE(mappos == snp0->mappos, file, lineno);

    for(i = 1; i < nSNPs; ++i) {
        snp = FileIndex_atNdx(fndx, i);
        SNPLoc_sanityCheck(snp, __FILE__, __LINE__);

        seekpos = FileIndex_getSeekPos(fndx, i);
        REQUIRE(seekpos >= snp0->seekpos, file, lineno);
        REQUIRE(seekpos == snp->seekpos, file, lineno);

        mappos = FileIndex_getMapPos(fndx, i);
        REQUIRE(mappos >= snp0->mappos, file, lineno);
        REQUIRE(mappos == snp->mappos, file, lineno);
        REQUIRE(i == FileIndex_atMapPos(fndx, mappos), file, lineno);

        snp0 = snp;
    }
}

#pragma GCC diagnostic pop

/** Construct an empty FileIndex */
FileIndex  *FileIndex_new(void) {
    FileIndex  *fndx = malloc(sizeof(FileIndex));

    checkmem(fndx, __FILE__, __LINE__);
    int         i;

    fndx->chunksAllocated = 0;
    fndx->nSNPs = 0;

    fndx->ptrsAllocated = NPTRS;
    fndx->snploc = malloc(fndx->ptrsAllocated * sizeof(SNPLoc));
    checkmem(fndx->snploc, __FILE__, __LINE__);
    for(i = 0; i < fndx->ptrsAllocated; ++i)
        fndx->snploc[i] = NULL;

    FileIndex_sanityCheck(fndx, __FILE__, __LINE__);
    fndx->nGtype = fndx->ploidy = 0;
    return fndx;
}

/** Free a FileIndex object */
void FileIndex_free(FileIndex * fndx) {
    long        i;

    if(fndx == NULL)
        return;
    for(i = 0; i < fndx->ptrsAllocated; ++i) {
        if(fndx->snploc[i] == NULL)
            break;
        free(fndx->snploc[i]);
        fndx->snploc[i] = NULL;
    }
    free(fndx->snploc);
    fndx->snploc = NULL;
    free(fndx);
    return;
}

/**
 * Reallocate the list of pointers so that the number of them is at
 * least as large as ptrsNeeded. All newly-allocated pointers are set
 * equal to NULL.
 */
void FileIndex_realloc(FileIndex * fndx, long ptrsNeeded) {
    myassert(fndx);
    long        i, goal;
    ldiv_t      quot;

    if(ptrsNeeded <= fndx->ptrsAllocated)
        return;

    /*
     * We allocate pointers in blocks of NPTRS pointers.
     * The first step is to figure out how many of pointers we wish to
     * have allocated after the call to realloc. This number is
     * "goal".
     */
    quot = ldiv(ptrsNeeded, (long) NPTRS);
    if(quot.rem != 0)
        goal = (quot.quot + 1) * NPTRS;
    else
        goal = quot.quot * NPTRS;

    /* realloc without the famous memory leak */
    SNPLoc    **tmp = realloc(fndx->snploc, goal * sizeof(fndx->snploc[0]));
    if(tmp == NULL)
        die("bad realloc", __FILE__, __LINE__);
    else
        fndx->snploc = tmp;

    for(i = fndx->ptrsAllocated; i < goal; ++i)
        fndx->snploc[i] = NULL;
    fndx->ptrsAllocated = goal;
    myassert(fndx->ptrsAllocated >= ptrsNeeded);
    FileIndex_sanityCheck(fndx, __FILE__, __LINE__);
    return;
}

/**
 * Record location data for a SNP in the FileIndex.
 *
 * @param[in] fndx points to FileIndex.
 * @param[in] seekpos seek position of new entry
 * @param[in] mappos map position (in cM) of new entry
 * @returns index of newly recorded SNP.
 */
long FileIndex_push(FileIndex * fndx, long seekpos, double mappos) {
    myassert(fndx);
    long        chunk, offset, i;

    chunk = fndx->nSNPs / CHUNKSIZE;
    offset = fndx->nSNPs % CHUNKSIZE;

    /* if necessary, use realloc to enlarge the array of pointers */
    if(chunk >= fndx->ptrsAllocated)
        FileIndex_realloc(fndx, chunk + 1);

    /* if necessary, allocate chunks of memory */
    while(chunk >= fndx->chunksAllocated) {
        i = fndx->chunksAllocated;
        myassert(i < fndx->ptrsAllocated);
        fndx->snploc[i] = malloc(CHUNKSIZE * sizeof(SNPLoc));
        checkmem(fndx->snploc[i], __FILE__, __LINE__);
        fndx->chunksAllocated += 1;
    }

    /* set data values for current SNP */
    fndx->snploc[chunk][offset].seekpos = seekpos;
    fndx->snploc[chunk][offset].mappos = mappos;

    ++fndx->nSNPs;

    return fndx->nSNPs;
}

/** Return the SNPLoc object for SNP with given index. **/
SNPLoc     *FileIndex_atNdx(FileIndex * fndx, long ndx) {
    ldiv_t      qr;

    myassert(fndx);
#ifndef NDEBUG
    if(ndx >= fndx->nSNPs) {
        fprintf(stderr,
                "ERR@%s:%d: ndx=%ld. Must not exceed fndx->nSNPs=%ld\n",
                __FILE__, __LINE__, ndx, fndx->nSNPs);
        dostacktrace(__FILE__, __LINE__, stderr);
        exit(EXIT_FAILURE);
    }
#endif
    qr = ldiv(ndx, (long) CHUNKSIZE);

    assert(qr.quot == ndx / CHUNKSIZE);
    assert(qr.rem == ndx % CHUNKSIZE);

    /* qr.quot is the index of the chunk */
    /* qr.rem is the index into that chunk */

    myassert(fndx->snploc[qr.quot]);

    return fndx->snploc[qr.quot] + qr.rem;
}

/** Return the map position (in cM) of SNP with given index. **/
double FileIndex_getMapPos(FileIndex * fndx, long ndx) {
    myassert(fndx);
    myassert(ndx < fndx->nSNPs);

    return FileIndex_atNdx(fndx, ndx)->mappos;
}

/** Return the seek position of SNP with given index. **/
long FileIndex_getSeekPos(FileIndex * fndx, long ndx) {
    myassert(fndx);
    myassert(ndx < fndx->nSNPs);

    return FileIndex_atNdx(fndx, ndx)->seekpos;
}

/** Return the number of SNPs. */
long FileIndex_nSNPs(FileIndex * fndx) {
    myassert(fndx);
    return fndx->nSNPs;
}

/**
 * Return index of smallest element in sorted array that is >=
 * mappos.
 */
long FileIndex_atMapPos(FileIndex * fndx, double mappos) {
    long        lo, mid, hi;
    double      hival, loval, midval;

#ifndef NDEBUG
    double      hm1val;
#endif

    lo = 0;
    hi = FileIndex_nSNPs(fndx) - 1;
    hival = FileIndex_atNdx(fndx, hi)->mappos;
    loval = FileIndex_atNdx(fndx, lo)->mappos;

    if(hival < loval) {
        fprintf(stderr, "WARNING@%s:%d hival (=%lg) < loval (=%lg)\n",
                __FILE__, __LINE__, hival, loval);
        return -1;
    }

    if(mappos > hival) {
        fprintf(stderr, "ERR@%s:%d mappos (=%lg) > hival (=%lg)\n",
                __FILE__, __LINE__, mappos, hival);
        return -1;
    }

    if(hi == 0)
        return 0;

    while(lo < hi) {
        mid = lo + (hi - lo) / 2;
        if(mid == lo)
            break;
        midval = FileIndex_atNdx(fndx, mid)->mappos;
        if(midval < mappos)
            lo = mid;
        else
            hi = mid;
    }
    loval = FileIndex_atNdx(fndx, lo)->mappos;
    if(loval >= mappos)
        hi = lo;

    hival = FileIndex_atNdx(fndx, hi)->mappos;

#ifndef NDEBUG
    if(hi > 0) {
        hm1val = FileIndex_atNdx(fndx, hi - 1)->mappos;
        myassert(hm1val < mappos);
    }
#endif

    myassert(hi < FileIndex_nSNPs(fndx));
    myassert(hival >= mappos);
    myassert(hi >= 0);

    return hi;
}

void FileIndex_printSummary(FileIndex * fndx, FILE * ofp) {
    fprintf(ofp, "FileIndex\n");
    fprintf(ofp, "  chunksAllocated=%ld ptrsAllocated=%ld nSNPs=%ld\n",
            fndx->chunksAllocated, fndx->ptrsAllocated, fndx->nSNPs);
    SNPLoc     *sl;

    if(fndx->nSNPs > 0) {
        sl = FileIndex_atNdx(fndx, 0);
        fprintf(ofp, "  SNPLoc[%5d]: seekpos=%ld mappos=%lg\n",
                0, sl->seekpos, sl->mappos);
    }
    if(fndx->nSNPs > 1) {
        long        last = fndx->nSNPs - 1;

        sl = FileIndex_atNdx(fndx, last);
        fprintf(ofp, "  SNPLoc[%5ld]: seekpos=%ld mappos=%lg\n",
                last, sl->seekpos, sl->mappos);
    }
    fprintf(ofp, "  ploidy=%u nGtype=%u\n", fndx->ploidy, fndx->nGtype);
}

void FileIndex_printChunks(FileIndex * fndx, FILE * ofp) {
    long        i, j;

    for(i = 0; i < fndx->ptrsAllocated; ++i) {
        if(fndx->snploc[i] != NULL) {
            fprintf(ofp, "Chunk %ld:\n", i);
            fprintf(ofp, " %10s %10s\n", "seekpos", "mappos");
            for(j = 0; j < CHUNKSIZE; ++j)
                fprintf(ofp, " %10ld %10lg\n",
                        fndx->snploc[i][j].seekpos,
                        fndx->snploc[i][j].mappos);
        }
    }
}

void FileIndex_print(FileIndex * fndx, FILE * ofp) {
    long        i;

    fprintf(ofp, "%10s %10s %10s\n", "ndx", "seekpos", "mappos");
    for(i = 0; i < FileIndex_nSNPs(fndx); ++i) {
        fprintf(ofp, "%10ld %10ld %10lg\n",
                i,
                FileIndex_getSeekPos(fndx, i), FileIndex_getMapPos(fndx, i));
    }
}

unsigned FileIndex_ploidy(const FileIndex *fndx) {
    return fndx->ploidy;
}

unsigned FileIndex_nGtype(const FileIndex *fndx) {
    return fndx->nGtype;
}

/**
 * Read through a file, storing the location of each SNP in a
 * structure of type FileIndex, which is then returned.
 *
 * @param[in] ifp Pointer to input stream.
 * @returns NULL if FILE pointer is NULL. If no sites are found, it
 * returns an empty FileIndex. Otherwise, function returns a pointer
 * to an object of type FileIndex.
 */
FileIndex  *FileIndex_readFile(FILE * ifp) {
    FileIndex  *fndx;
    int         rval = 0;
    double      mappos0 = -1.0;
    unsigned char sitedat[1000];
    int         ploidy = 1;

    if(ifp == NULL)
        return NULL;

    fndx = FileIndex_new();

    /*
     * Call Gtp_readHdr in order to position ifp at the start of the
     * data value. Get ploidy value from header.
     */
    rewind(ifp);
    Assignment *a = Gtp_readHdr(ifp);

    Assignment_setInt(a, "ploidy", &ploidy, MANDATORY);
    if(ploidy != 1 && ploidy != 2)
        eprintf("ERR@%s:%d: Bad ploidy in gtp file: %d\n",
                __FILE__,__LINE__,ploidy);
    Assignment_free(a);

    fndx->ploidy = ploidy;

    while(rval != EOF) {
        long        seekpos;
        double      mappos;

        seekpos = ftell(ifp);

        /*
         * The "NULL, 0" argument pairs cause Gtp_readSNP not to
         * return the snpId or the list of alleles.
         */
        rval = Gtp_readSNP(ifp, NULL, 0,    // don't get snpId
                           &mappos, NULL, 0,    // don't get alleles
                           sitedat, sizeof(sitedat), ploidy == 2);

        // Record the number of genotypes per SNP
        if(fndx->nGtype == 0 && rval > 0)
            fndx->nGtype = rval;

        if(mappos < mappos0) {
            fprintf(stderr,
                    "ERR@%s:%d: map positions are unsorted in file",
                    __FILE__, __LINE__);
            fprintf(stderr, " Detected at map position %14.12lg\n", mappos);
            exit(1);
        }
        if(Dbl_near(mappos, mappos0)) {
            /* ignore SNPs with duplicate mappos values */
            continue;
        }
        mappos0 = mappos;
        if(rval > 0)
            FileIndex_push(fndx, seekpos, mappos);
    }

    FileIndex_sanityCheck(fndx, __FILE__, __LINE__);
    return fndx;
}

/* Check sanity of a single ThreadBounds structure */
void ThreadBounds_sanityCheck(const ThreadBounds * tb, const char *file,
                              int lineno) {
    REQUIRE(tb != NULL, file, lineno);
    if(!(tb->seekpos[0] >= 0)) {
        ThreadBounds_print(tb, 1, stdout);
        fflush(stdout);
    }
    REQUIRE(tb->seekpos[0] >= 0, file, lineno);
    REQUIRE(tb->seekpos[1] >= tb->seekpos[0], file, lineno);
    REQUIRE(tb->seekpos[2] >= tb->seekpos[1], file, lineno);
    REQUIRE(tb->ndx[0] >= 0, file, lineno);
    REQUIRE(tb->ndx[1] >= tb->ndx[0], file, lineno);
    REQUIRE(tb->ndx[2] >= tb->ndx[1], file, lineno);
    REQUIRE(tb->mappos[0] >= 0, file, lineno);
    REQUIRE(tb->mappos[1] >= tb->mappos[0], file, lineno);
    REQUIRE(tb->mappos[2] >= tb->mappos[1], file, lineno);
    if(tb->ndx[1] > tb->ndx[0]) {
        REQUIRE(tb->seekpos[1] > tb->seekpos[0], file, lineno);
        REQUIRE(tb->mappos[1] > tb->mappos[0], file, lineno);
    }
    if(tb->seekpos[1] > tb->seekpos[0]) {
        REQUIRE(tb->ndx[1] > tb->ndx[0], file, lineno);
        REQUIRE(tb->mappos[1] > tb->mappos[0], file, lineno);
    }
    if(tb->mappos[1] > tb->mappos[0]) {
        if(!(tb->ndx[1] > tb->ndx[0])) {
            ThreadBounds_print(tb, 1, stdout);
            fflush(stdout);
        }
        REQUIRE(tb->ndx[1] > tb->ndx[0], file, lineno);
        REQUIRE(tb->seekpos[1] > tb->seekpos[0], file, lineno);
    }

    if(tb->ndx[2] > tb->ndx[1]) {
        REQUIRE(tb->seekpos[2] > tb->seekpos[1], file, lineno);
        REQUIRE(tb->mappos[2] > tb->mappos[1], file, lineno);
    }
    if(tb->seekpos[2] > tb->seekpos[1]) {
        REQUIRE(tb->ndx[2] > tb->ndx[1], file, lineno);
        REQUIRE(tb->mappos[2] > tb->mappos[1], file, lineno);
    }
    if(tb->mappos[2] > tb->mappos[1]) {
        REQUIRE(tb->ndx[2] > tb->ndx[1], file, lineno);
        REQUIRE(tb->seekpos[2] > tb->seekpos[1], file, lineno);
    }
}

/* Do sanity check on ThreadBounds structure of each thread. */
void ThreadBounds_arraySanityCheck(const ThreadBounds * tb, int nthreads,
                                   const char *file, int lineno) {
    const ThreadBounds *p;

    REQUIRE(tb != NULL, file, lineno);
    REQUIRE(nthreads > 0, file, lineno);
    for(p = tb; p < tb + nthreads; ++p)
        ThreadBounds_sanityCheck(p, file, lineno);
}

void ThreadBounds_print(const ThreadBounds * tb, int nthreads, FILE * ofp) {
    int         i;
    char        b1[100], b2[100], b3[100];

    fprintf(ofp, "#Bounds of Each Thread\n");
    fprintf(ofp, "#[initial, 1stFocal, lastFocal]\n");
    fprintf(ofp, "#[%26s][%26s][%26s]\n",
            strcenter("ndx", 26, b1, sizeof(b1)),
            strcenter("seekpos", 26, b2, sizeof(b2)),
            strcenter("mappos", 26, b3, sizeof(b3))
        );
    for(i = 0; i < nthreads; ++i) {
        fprintf(ofp, "#[%8ld,%8ld,%8ld]"
                "[%8ld,%8ld,%8ld]"
                "[%0.8lg,%0.8lg,%0.8lg]\n",
                tb[i].ndx[0],
                tb[i].ndx[1],
                tb[i].ndx[2],
                tb[i].seekpos[0],
                tb[i].seekpos[1],
                tb[i].seekpos[2],
                tb[i].mappos[0], tb[i].mappos[1], tb[i].mappos[2]);
    }
}

/**
 * Divide the genome into contiguous chunks of roughly equal size, so
 * that each chunk can be processed by a separate thread.
 *
 * On input, nthreads gives the number of threads, windowsize (the
 * number of SNPs in the sliding window), and ifp points to the input
 * file stream. The function returns a pointer to a newly-allocated
 * array of type ThreadBounds, which has an entry for each thread,
 * each of which is an object of type ThreadBounds. The i'th entry
 * defines the starting position for thread i.
 */
ThreadBounds *ThreadBounds_new(int nthreads,
                               double windowcm, FileIndex * fndx) {

    int         i;
    long        nSNPs = FileIndex_nSNPs(fndx);
    ThreadBounds *tb;

    tb = (ThreadBounds *) malloc(nthreads * sizeof(ThreadBounds));
    checkmem(tb, __FILE__, __LINE__);

    double      loCm = FileIndex_getMapPos(fndx, 0);

    /* perthread is # of focal SNPs per thread, rounded */
    long        perthread = floor(0.5 + nSNPs / ((double) nthreads));

    /*
     * These arrays are defined as in the documentation of the
     * ThreadBounds structure.
     */
    long        ndx[3], seekpos[3];
    double      mappos[3];

    for(i = 0; i < nthreads; ++i) {
        /* 1st focal SNP */
        ndx[1] = i * perthread;
        mappos[1] = FileIndex_getMapPos(fndx, ndx[1]);
        seekpos[1] = FileIndex_getSeekPos(fndx, ndx[1]);

        assert(ndx[1] < nSNPs);

        /* last focal SNP */
        ndx[2] = ndx[1] + perthread - 1;
        if(ndx[2] >= nSNPs)
            ndx[2] = nSNPs - 1;
        mappos[2] = FileIndex_getMapPos(fndx, ndx[2]);
        seekpos[2] = FileIndex_getSeekPos(fndx, ndx[2]);

        /* initial SNP: 1st to be compared with focal SNP */
        mappos[0] = mappos[1] - windowcm;
        if(mappos[0] < loCm)
            mappos[0] = loCm;
        ndx[0] = FileIndex_atMapPos(fndx, mappos[0]);
        if(ndx[0] < 0)
            ndx[0] = 0;
        mappos[0] = FileIndex_getMapPos(fndx, ndx[0]);
        seekpos[0] = FileIndex_getSeekPos(fndx, ndx[0]);

        assert(seekpos[0] >= 0);
        assert(mappos[0] >= 0);
        assert(ndx[0]==ndx[1] || windowcm>0.0);
        assert(mappos[0]==mappos[1] || windowcm>0.0);
        assert(seekpos[0]==seekpos[1] || windowcm>0.0);

        memcpy(tb[i].ndx, ndx, sizeof(ndx));
        memcpy(tb[i].seekpos, seekpos, sizeof(seekpos));
        memcpy(tb[i].mappos, mappos, sizeof(mappos));

#ifndef NDEBUG
        ThreadBounds_sanityCheck(tb + i, __FILE__, __LINE__);
#endif
    }
    if(tb[nthreads - 1].ndx[2] != nSNPs - 1) {
        tb[nthreads - 1].ndx[2] = nSNPs - 1;
        tb[nthreads - 1].mappos[2] = FileIndex_getMapPos(fndx, nSNPs - 1);
        tb[nthreads - 1].seekpos[2] = FileIndex_getSeekPos(fndx, nSNPs - 1);
    }
#ifndef NDEBUG
    ThreadBounds_arraySanityCheck(tb, nthreads, __FILE__, __LINE__);
#endif

    return tb;
}

void ThreadBounds_free(ThreadBounds * tb) {
    free(tb);
}

/** Return the index of the starting SNP. */
long ThreadBounds_ndx_initial(const ThreadBounds * tb) {
    return tb->ndx[0];
}

/** Return the index of the 1st focal SNP. */
long ThreadBounds_ndx_1stFocal(const ThreadBounds * tb) {
    return tb->ndx[1];
}

/** Return the index of the last focal SNP. */
long ThreadBounds_ndx_lastFocal(const ThreadBounds * tb) {
    return tb->ndx[2];
}

/** Return the seek position of the starting SNP. */
long ThreadBounds_seekpos_initial(const ThreadBounds * tb) {
    return tb->seekpos[0];
}

/** Return the seek position of the 1st focal SNP. */
long ThreadBounds_seekpos_1stFocal(const ThreadBounds * tb) {
    return tb->seekpos[1];
}

/** Return the seek position of the last focal SNP. */
long ThreadBounds_seekpos_lastFocal(const ThreadBounds * tb) {
    return tb->seekpos[2];
}

/** Return the map position of the starting SNP. */
double ThreadBounds_mappos_initial(const ThreadBounds * tb) {
    return tb->mappos[0];
}

/** Return the map position of the 1st focal SNP. */
double ThreadBounds_mappos_1stFocal(const ThreadBounds * tb) {
    return tb->mappos[1];
}

/** Return the map position of the last focal SNP. */
double ThreadBounds_mappos_lastFocal(const ThreadBounds * tb) {
    return tb->mappos[2];
}

