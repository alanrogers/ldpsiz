/**
@file preld.c
@anchor preld
@brief This file contains code related preld, a program that
predicts \f$\sigma_d^2\f$ from assumptions about population 
history.

`preld`, a program that predicts \f$\sigma_d^2\f$ from population history
==========================================================================

Parameter values, including those describing population history, may
be read either from the initialization file `ldpsiz.ini` or specified
on the command line. Command-line arguments override values in the
initialization file. Several methods are implemented for predicting
LD, and one or more of these may be specified either in the
initialization file or via the `--methods` argument. For the current list
of available methods, see \ref Model_alloc "Model_alloc".

Usage
-----

    usage: preld [options]
       where options may include:
       -u <x> or --mutation <x>
          set mutation rate/generation
       -b <x> or --nbins <x>
          specify the number of recombination rates
       -r <x> or --lo_r <x>
          low end of range of recombination rates in centimorgans
       -R <x> or --hi_r <x>
          high end of range of recombination rates in centimorgans
       -e or --equilibria
          show equilibrium for each epoch
       --log
          print log10 of sigdsq values
       -n or --twoNsmp
          set haploid sample size
       -E or --nextepoch
          move to next earlier epoch
       --twoN <x>
          set haploid pop size to x in current epoch
       -T <x> or --time <x>
          set length of current epoch to x generations
       -S or --printState
          print state vectors
       -m <method list> or --methods <method list>
          specify methods
       --exact
          don't use ODE to approximate difference equations
       -h or --help
          print this message

@copyright Copyright (c) 2014, Alan R. Rogers 
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include "misc.h"
#include "pophist.h"
#include "tokenizer.h"
#include "hill.h"
#include "hayes.h"
#include "strobeck.h"
#include "model.h"
#include "ini.h"

void        usage(void);
void        check_tn(double t, double n);

void usage(void) {
    fprintf(stderr, "usage: preld [options]\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-u <x> or --mutation <x>", "set mutation rate/generation");
    tellopt("-b <x> or --nbins <x>",
            "specify the number of recombination rates");
    tellopt("-r <x> or --lo_r <x>",
            "low end of range of recombination rates in centimorgans");
    tellopt("-R <x> or --hi_r <x>",
            "high end of range of recombination rates in centimorgans");
    tellopt("-e or --equilibria", "show equilibrium for each epoch");
    tellopt("--log", "print log10 of sigdsq values");
    tellopt("-n or --twoNsmp", "set haploid sample size");
    tellopt("-E or --nextepoch", "move to next earlier epoch");
    tellopt("--twoN <x>", "set haploid pop size to x in current epoch");
    tellopt("-T <x> or --time <x>",
            "set length of current epoch to x generations");
    tellopt("-S or --printState", "print state vectors");
    tellopt("-m <method list> or --methods <method list>", "specify methods");
    tellopt("--exact" , "don't use ODE to approximate difference equations");
    tellopt("-h or --help", "print this message");

    putc('\n', stderr);
    exit(1);
}

/**
 * Make sure values of t and N are valid.
 */
void check_tn(double t, double n) {
    if(isnan(t)) {
        fprintf(stderr, "\n\nERROR:");
        fprintf(stderr,
                "Specify --time <x> before --nextepoch on command line\n\n");
        usage();
    }
    if(isnan(n)) {
        fprintf(stderr, "\n\nERROR:");
        fprintf(stderr,
                "Specify --twoN <x> before --nextepoch on command line\n\n");
        usage();
    }
    return;
}

int main(int argc, char **argv) {

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"nbins", required_argument, 0, 'b'},
        {"exact", no_argument, 0, 'x'},
        {"log", no_argument, 0, 'l'},
        {"lo_r", required_argument, 0, 'r'},
        {"hi_r", required_argument, 0, 'R'},
        {"methods", required_argument, 0, 'm'},
        {"equilibria", no_argument, 0, 'e'},
        {"mutation", required_argument, 0, 'u'},
        {"printState", no_argument, 0, 'S'},
        {"twoNsmp", required_argument, 0, 'n'},
        {"nextepoch", no_argument, 0, 'E'},
        {"time", required_argument, 0, 'T'},
        {"twoN", required_argument, 0, 'N'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    /* parameters */
    int         nbins = 25;     /* number of recombination rates */
    double      lo_c = -1, hi_c = 0.003;
    int         doEquilibria = 0;
    int         doExact = 0;
    int         doLog = 0;
    int         printState = 0;
    double      u = 1e-4;
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;

    /* variables */
    int         twoNsmp = 0;
    ModelList  *ml = NULL;
    time_t      currtime = time(NULL);
    int         optndx, i, j, k, curr_epoch = 0, epochPending =
        0, phSetFromFile = 0;
    double      curr_t = strtod("NAN", 0), curr_twoN = strtod("NAN", 0);
    double      c, step_c, sigdsq;
    EpochLink  *linkedList = NULL;
    PopHist    *ph = NULL;
    const char *ofname = NULL;
    FILE       *ofp;
    char        buff[100] = { 0 }, method_str[100] = {
    0};

    printf("##########################################\n");
    printf("# preld: predict linkage disequilibrium #\n");
    printf("##########################################\n");

    putchar('\n');
#ifdef __TIMESTAMP__
    printf("# Program was compiled: %s\n", __TIMESTAMP__);
#endif
    printf("# Program was run     : %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    /* import definitions from initialization file */
    Ini        *ini = Ini_new(INIFILE);

    if(ini) {
        Ini_setInt(ini, "nbins", &nbins, !MANDATORY);
        if(Ini_setDbl(ini, "loCm", &lo_c, !MANDATORY))
            lo_c *= 0.01;
        if(Ini_setDbl(ini, "windowCm", &hi_c, !MANDATORY))
            hi_c *= 0.01;
        Ini_setInt(ini, "twoNsmp", &twoNsmp, !MANDATORY);
        if(Ini_setString
           (ini, "methods", method_str, sizeof(method_str), !MANDATORY))
            Ini_setInt(ini, "doEquilibria", &doEquilibria, !MANDATORY);
        Ini_setDbl(ini, "mutation", &u, !MANDATORY);
        if(Ini_setEpochLink(ini, &linkedList, !MANDATORY))
            phSetFromFile = 1;
        Ini_free(ini);
    }

    /* command line arguments */
    for(;;) {
        i = getopt_long(argc, argv, "b:r:R:m:n:eu:ET:Sh", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'b':
            nbins = strtod(optarg, 0);
            break;
        case 'r':
            /* multiply by 0.01 to convert cM to recombination rate */
            lo_c = 0.01 * strtod(optarg, 0);
            break;
        case 'R':
            /* multiply by 0.01 to convert cM to recombination rate */
            hi_c = 0.01 * strtod(optarg, 0);
            break;
        case 'm':
            /* parse comma-separated list of methods */
            snprintf(method_str, sizeof(method_str), "%s", optarg);
            break;
        case 'e':
            doEquilibria = 1;
            break;
        case 'x':
            doExact = 1;
            break;
        case 'l':
            doLog = 1;
            break;
        case 'u':
            myassert(optarg);
            u = strtod(optarg, 0);
            break;
        case 'n':
            myassert(optarg);
            twoNsmp = strtol(optarg, NULL, 10);
            break;
        case 'E':
            if(!epochPending) {
                fprintf(stderr, "Arg out of place: %s\n", "--nextEpoch");
                usage();
            }
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                phSetFromFile = 0;
                linkedList = NULL;
            }
            check_tn(curr_t, curr_twoN);
            linkedList = EpochLink_new(linkedList, curr_t, curr_twoN);
            ++curr_epoch;
            epochPending = 0;
            break;
        case 'T':
            myassert(optarg);
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                phSetFromFile = 0;
                linkedList = NULL;
            }
            curr_t = strtod(optarg, 0);
            epochPending = 1;
            break;
        case 'N':
            myassert(optarg);
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                phSetFromFile = 0;
                linkedList = NULL;
            }
            curr_twoN = strtod(optarg, 0);
            epochPending = 1;
            break;
        case 'S':
            printState = 1;
            break;
        case 'h':
        default:
            usage();
        }
    }

    if(argc > optind) {
        fprintf(stderr, "\nERROR: Extraneous arguments");
        for(i = optind; i < argc; ++i) {
            fprintf(stderr, " %s", argv[i]);
        }
        fputs("\n\n", stderr);
        usage();
    }

    /* Hill is default Model */
    if(method_str[0] == '\0')
        snprintf(method_str, sizeof(method_str), "%s", "Hill");
    ml = ModelList_alloc(method_str, twoNsmp);
    unsigned    nModels = ModelList_size(ml);

    if(lo_c <= 0.0)
        lo_c = hi_c / (2.0 * nbins);

    if(hi_c <= lo_c)
        die("highest recombination rate must exceed lowest.",
            __FILE__, __LINE__);

    if(epochPending && !phSetFromFile) {
        linkedList = EpochLink_new(linkedList, curr_t, curr_twoN);
    }

    /* default PopHist */
    if(linkedList == NULL) {
        /* epoch 0: t=300, n=1e5 */
        linkedList = EpochLink_new(linkedList, 300.0, 1e5);

        /* epoch 1: t=Inf, n=1e2 */
        linkedList = EpochLink_new(linkedList, strtod("Inf", 0), 1e2);
    }

    ph = PopHist_fromEpochLink(linkedList);

    if(ofname) {
        ofp = fopen(ofname, "w");
        if(ofp == NULL)
            eprintf("ERR@%s:%d: can't open %s for output",
                    __FILE__, __LINE__, ofname);
    } else
        ofp = stdout;

    step_c = (hi_c - lo_c) / (nbins - 1);

    if(twoNsmp > 0)
        fprintf(ofp, "# %-36s = %d\n", "Haploid sample size", twoNsmp);
    else
        fprintf(ofp, "# %-36s = %s\n", "Haploid sample size", "infinity");
    fprintf(ofp, "# %-36s = %s\n", "doEquilibria",
            (doEquilibria ? "yes" : "no"));
    fprintf(ofp, "# %-36s = %s\n", "Method",
            (doExact ? "No ODE" : "ODE for methods Hill and/or Strobeck"));
    fprintf(ofp, "# %-36s = %lg\n", "mutation_rate", u);
    fprintf(ofp, "# %-36s = %lg to %lg\n", "centimorgans",
            lo_c * 100.0, hi_c * 100.0);
    fprintf(ofp, "# %-36s = %d\n", "nbins", nbins);
#ifdef HAYES_MUTATION_ADJUSTMENT
    printf("# %-36s = %s\n", "Hayes mutation adjustment: E[rsq]", "1/(4Nc+2)");
#else
    printf("# %-36s = %s\n", "no Hayes mutation adjustment: E[rsq]",
           "1/(4Nc+1)");
#endif


    fprintf(ofp, "# %-36s =", "Models");
    for(i = 0; i < ModelList_size(ml); ++i) {
        fprintf(ofp, " %s", Model_lbl(ModelList_model(ml, i)));
        if(i + 1 < ModelList_size(ml))
            putc(',', ofp);
    }
    fputs("\n\n", ofp);

    PopHist_print_comment(ph, "#", ofp);
    putc('\n', ofp);

    /* Method label centered in a field of with lblWidth */
    int         fwidth = 11;    /* width of one field */
    int         lblWidth = fwidth;

    if(doEquilibria)
        lblWidth += PopHist_nepoch(ph) * (1 + lblWidth);

    char        fmt0[10], fmt1[10], fmt2[10];

    snprintf(fmt0, sizeof(fmt0), "#%%%ds", fwidth - 1);
    snprintf(fmt1, sizeof(fmt1), " %%%ds", fwidth);

    /* Header */
    fprintf(ofp, "#%s", strcenter("", fwidth - 1, buff, sizeof(buff)));
    for(i = 0; i < nModels; ++i) {
        const Model *model = ModelList_model(ml, i);

        int         currWidth = lblWidth;

        if(printState)
            currWidth += Model_stateDim(model) * (1 + fwidth);

        putc(' ', ofp);
        fprintf(ofp, "%s", strcenter(Model_lbl(model),
                                     currWidth, buff, sizeof(buff)));
    }
    putc('\n', ofp);

    fprintf(ofp, fmt0, "cM");
    for(i = 0; i < nModels; ++i) {
        if(doLog)
            fprintf(ofp, fmt1, "log10_LD");
        else
            fprintf(ofp, fmt1, "LD");

        const Model *model = ModelList_model(ml, i);

        for(j = 0; printState && j < Model_stateDim(model); ++j) {
            fprintf(ofp, fmt1, Model_stateLbl(model, j));
        }

        for(j = 0; doEquilibria && j < PopHist_nepoch(ph); ++j) {
            if(doLog)
                snprintf(buff, sizeof(buff), "log10_eq%01d", j);
            else
                snprintf(buff, sizeof(buff), "eq%01d", j);
            fprintf(ofp, fmt1, buff);
        }
    }
    putc('\n', ofp);

    snprintf(fmt0, sizeof(fmt0), "%%%d.%df", fwidth, fwidth - 2);
    if(doLog) {
        snprintf(fmt1, sizeof(fmt1), " %%%d.%df", fwidth, fwidth - 4);
        snprintf(fmt2, sizeof(fmt2), " %%%d.%dg", fwidth, fwidth - 7);
    } else {
        snprintf(fmt1, sizeof(fmt1), " %%%d.%df", fwidth, fwidth - 2);
        snprintf(fmt2, sizeof(fmt2), " %%%d.%dg", fwidth, fwidth - 5);
    }

    fflush(ofp);
    for(i = 0; i < nbins; ++i) {
        c = lo_c + i * step_c;  /* recombination rate */

        /* multiply 100 to convert from raw recombination to cM */
        fprintf(ofp, fmt0, 100.0 * c);

        for(j = 0; j < nModels; ++j) {
            const Model *model = ModelList_model(ml, j);
            ODE        *ode = ODE_new(model, odeAbsTol, odeRelTol);

            if(doExact)
                sigdsq = Model_exactLD(model, ode, c, u, ph);
            else
                sigdsq = ODE_ld(ode, c, u, ph);

            if(doLog)
                sigdsq = log10(sigdsq);
            fprintf(ofp, fmt1, sigdsq);

            for(k = 0; printState && k < Model_stateDim(model); ++k) 
                fprintf(ofp, fmt2, ODE_stateVal(ode, k));

            for(k = 0; doEquilibria && k < PopHist_nepoch(ph); ++k) {
                sigdsq = ODE_ldEq(ode, c, u, ph, k);
                if(doLog)
                    sigdsq = log10(sigdsq);
                fprintf(ofp, fmt1, sigdsq);
            }
            ODE_free(ode);
        }
        putc('\n', ofp);
    }

    PopHist_free(ph);
    if(linkedList)
        EpochLink_free(linkedList);
    ModelList_free(ml);

    return 0;
}
