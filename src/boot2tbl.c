/**
 * @file boot2tbl.c
 * @author Alan R. Rogers
 * @brief Read a block bootstrap file, as produced by obsld, and
 * print it in tabular form.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <assert.h>
#include "misc.h"
#include "boot.h"

void        usage(void);

void usage(void) {
    fprintf(stderr, "usage: boot2tbl [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-h or --help", "print this message");
    exit(1);
}

int main(int argc, char **argv) {

    Boot       *boot = NULL;
    char       *ifname = NULL;
    FILE       *ifp = NULL;
    int         bin, nBins, rep, nReps, i, optndx;

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    printf("###################################################\n"
           "# boot2tbl: convert bootstrap values into a table #\n"
           "###################################################\n");

    /* command line arguments */
    for(;;) {
        i = getopt_long(argc, argv, "h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'h':
        default:
            usage();
        }
    }

    /* remaining option gives file name */
    switch (argc - optind) {
    case 0:
        fprintf(stderr, "Command line must specify input file\n");
        usage();
        break;
    case 1:
        ifname = strdup(argv[optind]);
        ifp = fopen(ifname, "r");
        if(ifp == NULL)
            eprintf("ERR@%s:%d: Couldn't open \"%s\" for input",
                    __FILE__, __LINE__, ifname);
        printf("# input file: \"%s\"\n", ifname);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(ifname != NULL);
    assert(ifp != NULL);

    boot = Boot_restore(ifp);

    nBins = Boot_nBins(boot);
    nReps = Boot_nReps(boot);

    printf(" %4s %10s %10s %13s %13s %13s\n",
           "bin", "numerator", "denom", "mean_kb", "nobs", "rsqSum");
    for(rep = 0; rep < nReps; ++rep) {
        printf("# replicate %d\n", rep);
        for(bin = 0; bin < nBins; ++bin) {
            long        nobs;
            double      numerator, denominator, sep_kb, rsqSum;

            nobs = Boot_rawCounts(boot, rep, bin, &numerator,
                                  &denominator, &rsqSum, &sep_kb);

            printf(" %4d %10.8lg %10.8lg %13.8lg %13ld %13.8lg\n", bin,
                   numerator, denominator, sep_kb, nobs, rsqSum);
        }
    }
    fclose(ifp);
    free(ifname);
    Boot_free(boot);

    return 0;
}
