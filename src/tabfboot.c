/**
 * @file tabfboot.c
 * @author Alan R. Rogers
 * @brief Tabulate bootstrap parameter estimates produced by eld.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "boot.h"
#include "misc.h"
#include "pophist.h"
#include "tokenizer.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define NAMESIZE 30

void usage(void);

void usage(void) {
    fprintf(stderr, "usage: tabfboot [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");

    tellopt("-n <x> or --numTimes <x>", "number of time values");
    tellopt("--maxTime <x>", "max time value");
    tellopt("-q <q1> <q2> ... or --quantiles <q1> <q2> ...",
            "specify quantiles");
    tellopt ("-L or --log", "report 2N values as log10");
    tellopt("-h or --help", "print this message");
    fprintf(stderr,"\ninput should be an .fboot file produced by sald\n");
    exit(1);
}

int main(int argc, char **argv) {
    char ifname[200] = {'\0'};
    FILE *ifp = NULL;
    int i, j, nPr=0;
    int logScale = 0;
    double *prvals=NULL;

    /* MaxT is the largest value of t we want to tabulate. */
    double maxT = 1000.0; /* maximum t value to tabulate */
    int nT = 30;          /* number of t values to tabulate */

    /* command line arguments */
    i = 1;
    while(i < argc) {
        if(isopt("-h", "--help", argv[i]))
            usage();
        else if(isopt("-n", "--numTimes", argv[i])) {
            if(++i >= argc)
                usage();
             nT = strtol(argv[i], NULL, 10);
        }else if(isopt(NULL, "--maxTime", argv[i])) {
            if(++i >= argc)
                usage();
             maxT = strtod(argv[i], NULL);
        }else if(isopt("-L", "--log", argv[i])){
                logScale = true;
        }else if(isopt("-q", "--quantiles", argv[i])) {
            if(++i >= argc)
                usage();
            char *endptr;
            double q;
            /* count quantiles */
            for(nPr=0; i+nPr<argc; ++nPr) {
                q = strtod(argv[i+nPr], &endptr);
                if(endptr == argv[i+nPr]) /* argv not a number */
                    break;
                if(*endptr != '\0') {
                    fprintf(stderr, "%s:%d: Argument \"%s\" isn't a floating-point number",
                            __FILE__,__LINE__,argv[i+nPr]);
                    exit(EXIT_FAILURE);
                }
                if(!(q>0.0 && q<1.0)) {
                    fprintf(stderr, "%s:%d: Argument \"%s\" is outside range (0,1)",
                            __FILE__,__LINE__,argv[i+nPr]);
                    exit(EXIT_FAILURE);
                }
            }
            if(nPr == 0)
                usage();
            prvals = malloc(nPr * sizeof(prvals[0]));
            if(prvals==NULL)
                die("malloc", __FILE__,__LINE__);
            for(j=0; j<nPr; ++j) {
                q = strtod(argv[i+j], &endptr);
                if(endptr == argv[i+j]) /* argv not a number */
                    die("This can't happen",__FILE__,__LINE__);
                prvals[j] = q;
            }
            i += nPr-1;
        }else if(argv[i][0] == '-') {
            fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
            usage();
        }else {
            if(ifp != NULL)
                eprintf("ERR@%s:%d: only one input file allowed",
                        __FILE__, __LINE__);
            snprintf(ifname, sizeof(ifname), "%s", argv[i]);
            ifp = fopen(ifname, "r");
            if(ifp == NULL)
                eprintf("ERR@%s:%d: Couldn't open %s for input",
                        __FILE__, __LINE__, ifname);
        }
        ++i;
    }
    
    if(ifp == NULL)
        usage();

    /* generate vector, tvals, of time values */
    double w;
    if(maxT > 0.0) {
        if(nT < 2)
            nT = 2;
        w = maxT/(nT-1);
    }else{
        nT = 1;
        maxT = 0.0;
        w = 1.0;
    }
    double *tvals = malloc(nT*sizeof(tvals[0]));
    checkmem(tvals, __FILE__, __LINE__);
    for(i=0; i < nT; ++i)
        tvals[i] = i*w;

    int nparams=0, ntokens;
    
    char tokbuff[500], buff[500];
    Tokenizer *tkz = Tokenizer_new(100);

    while(nparams==0) {
        if(NULL == fgets(tokbuff, (int) sizeof(tokbuff), ifp))
            die("Couldn't read data", __FILE__, __LINE__);

        if(!strchr(tokbuff, '\n'))
            eprintf("ERR@%s:%d: input tokbuffer overflow."
                    " tokbuff size: %d\n", __FILE__, __LINE__, sizeof(tokbuff));

        Tokenizer_split(tkz, tokbuff, " ");
        nparams = Tokenizer_strip(tkz, " \n");
    }

    /* 1st pass: count datasets only */
    int rep, nreps=0;
    if(ifp==NULL)
        die("ifp==NULL", __FILE__,__LINE__);
    for(rep=0;
        NULL != fgets(buff, (int) sizeof(buff), ifp);
        ++rep) {
        if(!strchr(buff, '\n'))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));
        ++nreps;
    }
    rewind(ifp);

    // Read header. Determine whether fboot file tabulates twoNinv or
    // twoN. If it tabulates twoNinv, then set invert2N=true. Otherwise,
    // set invert2N=false.
    int invert2N = false;

    if(NULL==fgets(tokbuff, (int) sizeof(tokbuff), ifp))
        die("Unexpected EOF", __FILE__, __LINE__); // read header
    Tokenizer_split(tkz, tokbuff, " ");
    ntokens = Tokenizer_strip(tkz, " \n");
    if(ntokens != nparams)
        eprintf("ERR@%s:%s:%d: ntokens=%d; should equal nparams=%d",
                __func__,__FILE__,__LINE__,ntokens,nparams);
    if(0 != strncmp("twoN", Tokenizer_token(tkz,0), 4))
        eprintf("%s:%s:%d: read %s; expecting %s or %s",
                __FILE__,__func__,__LINE__,
                Tokenizer_token(tkz,0), "twoN", "twoNinv");
    if(0 == strncmp("twoNinv", Tokenizer_token(tkz,0),7))
        invert2N = true;

    int nepoch = 1 + nparams/2;
    PopHist *ph = PopHist_newEmpty(nepoch);
    double **twoN = malloc(nT * sizeof(twoN[0]));
    checkmem(twoN, __FILE__, __LINE__);
    for(i=0; i<nT; ++i) {
        twoN[i] = malloc(nreps * sizeof(twoN[0][0]));
        checkmem(twoN[i], __FILE__, __LINE__);
        memset(twoN[i], 0, nreps*sizeof(twoN[0][0]));
    }

    /* 2nd pass: fill twoN array */
    for(rep=0;
        NULL != fgets(tokbuff, (int) sizeof(tokbuff), ifp);
        ++rep) {

        if(rep >= nreps)
            eprintf("%s:%d: rep=%d is equal to nreps=%d\n",
                    __FILE__,__LINE__,rep,nreps);
        if(!strchr(tokbuff, '\n'))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(tokbuff));

        /* Fill PopHist with estimates from current dataset */
        Tokenizer_split(tkz, tokbuff, " ");
        ntokens = Tokenizer_strip(tkz, " \n");
        if(ntokens != nparams)
            eprintf("ERR@%s:%s:%d: ntokens=%d; should equal nparams=%d",
                    __func__,__FILE__,__LINE__,ntokens,nparams);
        double duration, curr2N;
        for(i=0; i< nepoch-1; ++i) {
            if(invert2N)
                curr2N = 1.0/strtod(Tokenizer_token(tkz, 2*i), NULL);
            else
                curr2N = strtod(Tokenizer_token(tkz, 2*i), NULL);
            duration = strtod(Tokenizer_token(tkz, 2*i+1), NULL);

            PopHist_setTwoN(ph, i, curr2N);
            PopHist_setDuration(ph, i, duration);
        }
        if(invert2N)
            curr2N = 1.0/strtod(Tokenizer_token(tkz, 2*(nepoch-1)), NULL);
        else
            curr2N = strtod(Tokenizer_token(tkz, 2*(nepoch-1)), NULL);
        PopHist_setTwoN(ph, nepoch-1, curr2N);

        for(i=0; i<nT; ++i) {
            j = PopHist_findEpoch(ph, tvals[i]);
            curr2N = PopHist_twoN(ph, j);
            twoN[i][rep] = curr2N;
        }
    }

    if(nPr == 0) {
        nPr = 5;
        prvals = malloc(nPr * sizeof(prvals[0]));
        checkmem(prvals, __FILE__,__LINE__);
        prvals[0] = 0.025;
        prvals[1] = 0.25;
        prvals[2] = 0.5;
        prvals[3] = 0.75;
        prvals[4] = 0.975;
    }

    printf("%8s", "t");
    for(j=0; j<nPr; ++j) {
        snprintf(buff, sizeof(buff), "%sq%lg", (logScale ? "log" : ""), prvals[j]);
        printf(" %9s", buff);
    }
    snprintf(buff, sizeof(buff), "%s2Nmean", (logScale ? "log" : ""));
    printf(" %9s\n", buff);

    for(i=0; i < nT; ++i) {

        printf("%8.3lf", tvals[i]);

        qsort(twoN[i], nreps, sizeof(twoN[i][0]), compareDoubles);
        for(j=0; j<nPr; ++j) {
            double q = interpolate(prvals[j], &twoN[i][0], nreps);
            printf(" %9.5lf", (logScale ? log10(q) : q));
        }

        double twoNmean = 0.0;
        for(j=0; j<nreps; ++j)
            twoNmean += twoN[i][j];
        twoNmean /= nreps;
        printf(" %9.5lf", (logScale ? log10(twoNmean) : twoNmean));
        putchar('\n');
    }
    Tokenizer_free(tkz);
	PopHist_free(ph);
	fclose(ifp);
    for(i=0; i<nT; ++i)
        free(twoN[i]);
    free(twoN);
    free(tvals);
    return 0;
}
