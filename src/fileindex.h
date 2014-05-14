/**
 * @file fileindex.h
 * @author Alan R. Rogers
 * @brief Header for fileindex.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef ESTLD_FILEINDEX_H
#define ESTLD_FILEINDEX_H

#include <stdio.h>
#include "typedefs.h"

/**
 * Each ThreadBounds object contains the information needed by an
 * individual thread to figure out where in the input file to begin
 * reading SNPs and which SNP is its first focal SNP. The definition
 * is exposed here so that it will be possible to do array indexing on
 * an array of objects of type ThreadBounds.
 */
struct ThreadBounds {
    /*
     * A ThreadBounds object consists of 3 arrays, each with 3
     * entries. The arrays are
     *
     * ndx    : the (base-zero) indices of SNPs
     * seekpos: position in file, as reported by ftell
     * mappos : map position on chromosome, in centimorgans
     *
     * With each array, entry 0 refers to the the first SNP to be read
     * by this thread, entry 1 to the first focal SNP, and entry 2 to
     * the last focal SNP. 
     */
    long        ndx[3], seekpos[3];
    double      mappos[3];
};

FileIndex  *FileIndex_new(void);
void        FileIndex_free(FileIndex * fndx);
long        FileIndex_push(FileIndex * fndx, long seekpos, double mappos);
SNPLoc     *FileIndex_atNdx(FileIndex * fndx, long ndx);
long        FileIndex_atMapPos(FileIndex * fndx, double mappos);
long        FileIndex_nSNPs(FileIndex * fndx);
void        FileIndex_printSummary(FileIndex * fndx, FILE * ofp);
void        FileIndex_printChunks(FileIndex * fndx, FILE * ofp);
void        FileIndex_print(FileIndex * fndx, FILE * ofp);
FileIndex  *FileIndex_readFile(FILE * ifp);
void        FileIndex_realloc(FileIndex * fndx, long ptrsNeeded);
double      FileIndex_getMapPos(FileIndex * fndx, long ndx);
long        FileIndex_getSeekPos(FileIndex * fndx, long ndx);
void        FileIndex_sanityCheck(FileIndex * fi, const char *file,
                                  int lineno);
unsigned    FileIndex_nPtrs(void);
unsigned    FileIndex_chunkSize(void);
long        FileIndex_nPtrsAllocated(const FileIndex *fndx);

long        SNPLoc_seekpos(const SNPLoc *sl);
double      SNPLoc_mappos(const SNPLoc *sl);

ThreadBounds *ThreadBounds_new(int nthreads, double windowcm,
                               FileIndex * fndx);
void        ThreadBounds_free(ThreadBounds * tb);
void        ThreadBounds_print(const ThreadBounds * tb, int nthreads,
                               FILE * ofp);
void        ThreadBounds_sanityCheck(const ThreadBounds * tb,
                                     const char *file, int lineno);
void        ThreadBounds_arraySanityCheck(const ThreadBounds * tb,
                                          int nthreads, const char *file,
                                          int lineno);
long        ThreadBounds_ndx_initial(const ThreadBounds * tb);
long        ThreadBounds_ndx_1stFocal(const ThreadBounds * tb);
long        ThreadBounds_ndx_lastFocal(const ThreadBounds * tb);
long        ThreadBounds_seekpos_initial(const ThreadBounds * tb);
long        ThreadBounds_seekpos_1stFocal(const ThreadBounds * tb);
long        ThreadBounds_seekpos_lastFocal(const ThreadBounds * tb);
double      ThreadBounds_mappos_initial(const ThreadBounds * tb);
double      ThreadBounds_mappos_1stFocal(const ThreadBounds * tb);
double      ThreadBounds_mappos_lastFocal(const ThreadBounds * tb);

#endif
