/**
 * @file fithayes.c
 * @author Alan R. Rogers
 * @brief Functions implementing the method of Hayes et al (Genome
 * Research, 13(4): 635-643).
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <float.h>
#include <getopt.h>
#include <pthread.h>
#include <string.h>
#include "misc.h"
#include "tokenizer.h"
#include "assign.h"
#include "ini.h"

#define NTRIES 20
#define MINKB  1.0
#define MAXKB  300.0

void        usage(void);
void        write_data(FILE * outfile, int nvars, int ncc, double *cc,
                       double **ld, const char **lbl);
int         read_data(FILE * ifp, int nbins, double *cm, double *rsq,
                      long *nobs);

void usage(void) {
    fprintf(stderr, "usage: fithayes [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
#if 0
    tellopt("-n or --twoNsmp", "set haploid sample size");
#endif
    tellopt("-h or --help", "print this message");

    exit(1);
}

void
write_data(FILE * outfile, int nvars, int ncc, double *cc,
           double **ld, const char **lbl) {
    int         i, j, nlbl = nvars + 1;

    fprintf(outfile, "#nbins : %d\n", ncc);
    fprintf(outfile, "#%7s", lbl[0]);
    for(i = 1; i < nlbl; ++i) {
        fprintf(outfile, " %14s", lbl[i]);
    }
    putc('\n', outfile);
    for(i = 0; i < ncc; ++i) {
        fprintf(outfile, "%8.5f", 100.0 * cc[i]);
        for(j = 0; j < nvars; ++j) {
            fprintf(outfile, " %14.5f", ld[j][i]);
        }
        putc('\n', outfile);
    }
}

/**
 * Read the data file produced by eld.
 *
 * @param[in] ifp Points to file produced by eld.
 * @param[in] nbins The length of all arrays.
 * @param[out] cm An array giving the average separation (in centimorgans)
 * between pairs of SNPs within the various bins.
 * @param[out] rsq An array of estimates of r^2.
 * @param[out] nobs An array whose i'th entry is the number of pairs
 * of SNPs contributing to bin i.
 *
 * @returns number of lines read
 */
int read_data(FILE * ifp, int nbins, double *cm, double *rsq, long *nobs) {
    char        buff[200];
    int         inData = 0, ntokens, i, tokensExpected = 0;
    Tokenizer  *tkz = Tokenizer_new(50);

    /* skip until beginning of data */
    while(!inData && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        Tokenizer_split(tkz, buff, " \t");  /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens < 3)
            continue;
        if(strcmp(Tokenizer_token(tkz, 0), "cM") == 0) {
            inData = 1;
            switch (ntokens) {
            case 4:
                /* fall through */
            case 6:
                tokensExpected = ntokens;
                break;
            default:
                fprintf(stderr, "Current tokens:");
                Tokenizer_print(tkz, stderr);
                eprintf("ERR@%s:%d: got %d tokens rather than 4 or 6",
                        __FILE__, __LINE__, ntokens);
            }
        }
    }

    if(!inData)
        die("Couldn't find data in input file", __FILE__, __LINE__);

    /* read data */
    i = 0;
    while(i < nbins && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        ntokens = Tokenizer_split(tkz, buff, " \t");    /* tokenize */
        if(ntokens == 0)
            continue;

        if(ntokens != tokensExpected)
            eprintf("ERR@%s:%d: got %d tokens rather than %d",
                    __FILE__, __LINE__, ntokens, tokensExpected);

        cm[i] = strtod(Tokenizer_token(tkz, 0), NULL);
        nobs[i] = strtol(Tokenizer_token(tkz, 2), NULL, 10);
        rsq[i] = strtod(Tokenizer_token(tkz, tokensExpected-1), NULL);

        ++i;
    }

    Tokenizer_free(tkz);
    tkz = NULL;
    return i;                   /* return number of lines read */
}

int main(int argc, char **argv) {

    int         optndx;

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
#if 0
        {"twoNsmp", required_argument, 0, 'n'},
#endif
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    double     *cm, *rsq_obs;

    int         nbins;
    int         twoNsmp = 0;    /* number of haploid samples */

    time_t      currtime = time(NULL);
    char        method[20];
    char       *ifname = NULL;
    FILE       *ifp = NULL;
    int         i, rval;
    long int   *nobs;

    printf("##################################################\n"
           "# fithayes: fit population parameters to linkage #\n"
           "#           disequilibrium curve                 #\n"
           "##################################################\n");

    putchar('\n');
#ifdef __TIMESTAMP__
    printf("# Program was compiled: %s\n", __TIMESTAMP__);
#endif
    printf("# Program was run     : %s\n", ctime(&currtime));

    printf("# fithayes command:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

#if 0
    /* import definitions from initialization file */
    Ini        *ini = Ini_new(INIFILE);

    if(ini) {
        Ini_setInt(ini, "twoNsmp", &twoNsmp, !MANDATORY);
        Ini_free(ini);
        ini = NULL;
    }
#endif

    /* command line arguments */
    for(;;) {
        i = getopt_long(argc, argv, "h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
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
            eprintf("ERR@%s:%d: Couldn't open %s for input",
                    __FILE__, __LINE__, ifname);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(ifname != NULL);
    assert(ifp != NULL);

    /* specify model */
    snprintf(method, sizeof(method), "Hayes");

    /* read assignment statements in input file */
    Assignment *asmt = Assignment_readEld(ifp);

    Assignment_setInt(asmt, "nbins", &nbins, MANDATORY);
    Assignment_setInt(asmt, "Haploid sample size", &twoNsmp, MANDATORY);
    Assignment_free(asmt);
    asmt = NULL;

    cm = (double *) malloc(nbins * sizeof(cm[0]));
    checkmem(cm, __FILE__, __LINE__);

    rsq_obs = (double *) malloc(nbins * sizeof(rsq_obs[0]));
    checkmem(rsq_obs, __FILE__, __LINE__);

    nobs = (long int *) malloc(nbins * sizeof(nobs[0]));
    checkmem(nobs, __FILE__, __LINE__);

    rval = read_data(ifp, nbins, cm, rsq_obs, nobs);
    if(rval != nbins)
        eprintf("ERR@%s:%d: Couldn't read %d lines of data from \"%s\"."
                " Only found %d lines", nbins, ifname, rval);
    printf("# %-35s = %d\n", "Haploid sample size", twoNsmp);

#ifdef HAYES_MUTATION_ADJUSTMENT
        printf("Adjusting for mutation: E[rsq] = 1/(4Nc+2)\n");
#else
        printf("Not adjusting for mutation: E[rsq] = 1/(4Nc+1)\n");
#endif

    /*
     * Use Hayes model to estimate twoN for each bin.
     *
     * Assume sigdsq ~= rsq.
     * Subtract 1/twoNsmp to account for sampling bias.
     * Solve rsq = 1/(1 + 2*twoN*c):
     * twoN = ((1/rsq) -1)/(2*c)
     */
    printf("#%10s %9s\n", "generation", "twoN");
    for(i = nbins - 1; i >= 0; --i) {
        double      c = 0.01 * cm[i];   /* cM --> recombination rate */
        double      rsq = rsq_obs[i] - 1.0 / twoNsmp;

        if(rsq < 0.0)
            rsq = 0.0;
#ifdef HAYES_MUTATION_ADJUSTMENT
        /* accounts for mutation */
        double      twoN = ((1.0 / rsq) - 2.0) / (2.0 * c);
#else
        /* no mutation */
        double      twoN = ((1.0 / rsq) - 1.0) / (2.0 * c);
#endif
        double      t = 1.0 / (2.0 * c);

        printf("%11.1lf %10.1lf\n", t, twoN);
    }

    if(ifname)
        free(ifname);
    if(ifp)
        fclose(ifp);
    free(cm);
    free(rsq_obs);
    free(nobs);

    return 0;
}
