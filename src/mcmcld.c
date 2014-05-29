/**
@file mcmcld.c
@page mcmcld
@anchor mcmcld
@brief `mcmcld` uses Markov-Coupled Markov-Chain Monte Carlo (MCMCMC)
to estimate the cost function \f$\sigma_d^2\f$.

mcmcld, MCMCMC estimates of population history from \f$\sigma_d^2\f$
====================================================================

This program needs work. Don't trust the current version.

@copyright Copyright (c) 2014, Alan R. Rogers 
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <pthread.h>
#include <string.h>
#include "misc.h"
#include "boot.h"
#include "pophist.h"
#include "tokenizer.h"
#include "ini.h"
#include "model.h"
#include "assign.h"
#include "chain.h"

void        usage(void);
int         read_data(FILE * ifp, int nbins, int twoNsmp, double *cm,
                      double *sigdsq);
int         isopt(const char *shortOpt, const char *longOpt, const char *arg);

void usage(void) {
    fprintf(stderr, "usage: mkfitld [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-u <x> or --mutation <x>", "set mutation rate/generation");
    tellopt("--dNinv <x>", "width of proposal distribution for 1/2N");
    tellopt("--dt <x>", "width of proposal distribution for t");
    tellopt("-n or --twoNsmp", "set haploid sample size");
    tellopt("--twoN <x>", "set haploid pop size to x in current epoch");
    tellopt("-T <x> or --time <x>",
            "set length of current epoch to x generations");
    tellopt("-t <x> or --temp <x>",
            "set temperature difference between adjacent chains");
    tellopt("-E or --nextepoch", "move to next earlier epoch");
    tellopt("-c <x> or --chains <x>", "number of parallel chains");
    tellopt("--nreps <x>", "specify max iterations for Markov Chain");
    tellopt("-I <x> or --outInterval <x>", "iterations between output lines");
    tellopt("-S <x> or --swapInterval <x>",
            "iterations between attempts to swap adjacent chains");
    tellopt("-m <method> or --methods <method>",
            "specify method \"Hill\" or \"Strobeck\"");
    tellopt("-h or --help", "print this message");

    exit(1);
}

/**
 * Log of objective function: negative of sum of squared differences
 * between observed and expected sigma_d^2 vectors.
 *
 * @param[in] ph Current population history
 * @param[in] u mutation rate per site per generation
 * @param[in] nbins Number of values in vectors obs and c
 * @param[in] sigdsq Vector of nbins values, the observed values of
 * sigdsq.
 * @param[in] c Vector of nbins values, the recombination rates
 * associated with the values in sigdsq.
 */
double lnObjFun(PopHist * ph, double u, ODE * ode,
                int nbins, double *sigdsq, double *c) {
    double      badness;
    double      exp_sigdsq[nbins];
    int         i;

    /* get vector of expected values of sigdsq */
    ODE_ldVec(ode, exp_sigdsq, nbins, c, u, ph);

    badness = 0.0;
    for(i = 0; i < nbins; ++i) {
        double      diff = exp_sigdsq[i] - sigdsq[i];

        badness += diff * diff;
    }

    if(!isfinite(badness)) {
#ifdef DEBUG
        for(i = 0; i < nbins; ++i) {
            if(!isfinite(exp_sigdsq[i])) {
                fprintf(stderr, "WARNING@%s:%d: exp_sigdsq[%d]=%lg\n",
                        __FILE__, __LINE__, i, exp_sigdsq[i]);
            }
            if(!isfinite(sigdsq[i]))
                eprintf("ERR@%s:%d: obs[%d]=%lg",
                        __FILE__, __LINE__, i, sigdsq[i]);
        }
#endif
        badness = DBL_MAX;
    }

    return -badness;
}

/**
 * Read the data file produced by eld.
 *
 * @param[in] ifp Points to file produced by eld.
 * @param[in] nbins The length of all arrays.
 * @param[out] twoNsmp Number of gene copies sampled per locus.
 * @param[out] cm An array giving the average separation (in centimorgans)
 * between pairs of SNPs within the various bins.
 * @param[out] sigdsq An array of estimates of sigma_d^2.
 *
 * @returns number of lines read
 */
int read_data(FILE * ifp, int nbins, int twoNsmp, double *cm, double *sigdsq) {
    char        buff[200];
    int         inData = 0, ntokens, i, tokensExpected = 0;
    Tokenizer  *tkz = Tokenizer_new(50);

    rewind(ifp);
    /* skip until beginning of data */
    while(!inData && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        Tokenizer_split(tkz, buff, " \t");  /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens < 2)
            continue;
        if(strcmp(Tokenizer_token(tkz, 0), "cM") == 0) {
            inData = 1;
            switch (ntokens) {
            case 2:
                /* fall through */
            case 3:
                /* fall through */
            case 5:
                tokensExpected = ntokens;
                break;
            default:
                fprintf(stderr, "Current tokens:");
                Tokenizer_print(tkz, stderr);
                eprintf("ERR@%s:%d: got %d tokens rather than 2, 3, or 5",
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
        sigdsq[i] = strtod(Tokenizer_token(tkz, 1), NULL);

        ++i;
    }

    Tokenizer_free(tkz);
    tkz = NULL;
    return i;                   /* return number of lines read */
}

int main(int argc, char **argv) {

    double      u = 1e-6;
    double     *cc, *sigdsq_obs;

    int         nbins;
    int         twoNsmp;        /* number of haploid samples */
    int         nchains = 1;    /* number of parallel chains */
    unsigned    outInterval = 1;    /* iterations between output lines */
    unsigned    swapInterval = 1;   /* iterations between attempts to swap */
    int         curr_epoch = 0, epochPending = 1, phSetFromFile = 0;
    double      curr_t = strtod("Inf", 0);
    double      curr_twoN = 1000.0;
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;

    char        method[20] = { 0 };
    EpochLink  *linkedList = NULL;
    char       *ifname = NULL;
    FILE       *ifp = NULL;
    int         i, rval;
    unsigned    rep, nreps = 100;
    PopHist    *ph_init;
    time_t      currtime = time(NULL);

    double      dt = 5e2;       /* width of proposal distribution: t    */
    double      dNinv = 0.05;   /* width of proposal distribution: 1/2N */
    double      temp = 0.7;     /* temperature factors for parallel chains */

    printf("###############################################\n"
           "# mcmcld: Use MCMCMC to fit population        #\n"
           "# parameters to linkage disequilibrium curve  #\n"
           "###############################################\n");
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
        Ini_setString(ini, "methods", method, sizeof(method), !MANDATORY);
        Ini_setDbl(ini, "mutation", &u, !MANDATORY);
        Ini_setDbl(ini, "temperature", &temp, !MANDATORY);
        Ini_setUnsignedInt(ini, "output interval", &outInterval, !MANDATORY);
        Ini_setUnsignedInt(ini, "swap interval", &swapInterval, !MANDATORY);
        Ini_setUnsignedInt(ini, "MCMC iterations", &nreps, !MANDATORY);
        if(Ini_setEpochLink(ini, &linkedList, !MANDATORY))
            phSetFromFile = 1;
        Ini_free(ini);
        ini = NULL;
    }

    /* command line arguments */
    i = 1;
    while(i < argc) {
        if(isopt("-h", "--help", argv[i]))
            usage();
        else if(isopt("-c", "--chains", argv[i])) {
            if(++i >= argc)
                usage();
            nchains = strtol(argv[i], NULL, 10);
        } else if(isopt("-E", "--nextepoch", argv[i])) {
            if(!epochPending) {
                fprintf(stderr, "Arg out of place: %s\n", "--nextEpoch");
                usage();
            }
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                linkedList = NULL;
                phSetFromFile = 0;
            }
            linkedList = EpochLink_new(linkedList, curr_t, curr_twoN);
            ++curr_epoch;
            epochPending = 0;
        } else if(isopt("-I", "--outInterval", argv[i])) {
            if(++i >= argc)
                usage();
            outInterval = strtoul(argv[i], NULL, 10);
        } else if(isopt("-m", "--methods", argv[i])) {
            if(++i >= argc)
                usage();
            snprintf(method, sizeof(method), "%s", argv[i]);
        } else if(isopt("-N", "--twoN", argv[i])) {
            if(++i >= argc)
                usage();
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                linkedList = NULL;
                phSetFromFile = 0;
            }
            curr_twoN = strtod(argv[i], 0);
            epochPending = 1;
        } else if(isopt("-n", "--twoNsmp", argv[i])) {
            if(++i >= argc)
                usage();
            twoNsmp = strtol(argv[i], NULL, 10);
        } else if(isopt("-r", "--nreps", argv[i])) {
            if(++i >= argc)
                usage();
            nreps = strtol(argv[i], NULL, 10);
        } else if(isopt("-S", "--swapInterval", argv[i])) {
            if(++i >= argc)
                usage();
            swapInterval = strtoul(argv[i], NULL, 10);
        } else if(isopt("-T", "--time", argv[i])) {
            if(++i >= argc)
                usage();
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                linkedList = NULL;
                phSetFromFile = 0;
            }
            curr_t = strtod(argv[i], 0);
            epochPending = 1;
        } else if(isopt("-t", "--temp", argv[i])) {
            if(++i >= argc)
                usage();
            temp = strtod(argv[i], 0);
        } else if(isopt("-u", "--mutation", argv[i])) {
            if(++i >= argc)
                usage();
            u = strtod(argv[i], 0);
        } else if(isopt("-W", "--dNinv", argv[i])) {
            if(++i >= argc)
                usage();
            dNinv = strtod(argv[i], 0);
        } else if(isopt("-w", "--dt", argv[i])) {
            if(++i >= argc)
                usage();
            dt = strtod(argv[i], 0);
        } else if(argv[i][0] == '-') {
            fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
            usage();
        } else {
            if(ifp != NULL)
                eprintf("ERR@%s:%d: only one input file allowed",
                        __FILE__, __LINE__);
            ifname = strdup(argv[i]);
            ifp = fopen(ifname, "r");
            if(ifp == NULL)
                eprintf("ERR@%s:%d: Couldn't open %s for input",
                        __FILE__, __LINE__, ifname);
        }
        ++i;
    }
    if(ifp == NULL) {
        fprintf(stderr, "ERR@%s:%d: no input file", __FILE__, __LINE__);
        usage();
    }
    assert(ifname != NULL);
    assert(ifp != NULL);

    /* specify default model */
    if(method[0] == '\0')
        snprintf(method, sizeof(method), "Hill");

    /* If more than one method was specified, use only first. */
    char       *p = strchr(method, ',');

    if(p) {
        fprintf(stderr, "WARNING@%s:%d: Models=\"%s\".",
                __FILE__, __LINE__, method);
        *p = '\0';
        fprintf(stderr, " Using only \"%s\".\n", method);
    }

    if(epochPending && !phSetFromFile)
        linkedList = EpochLink_new(linkedList, curr_t, curr_twoN);

    ph_init = PopHist_fromEpochLink(linkedList);

    if(PopHist_nepoch(ph_init) < 1) {
        fprintf(stderr, "Initial PopHist parameters unspecified\n");
        usage();
    }

    printf("# Initial PopHist:\n");
    PopHist_print_comment(ph_init, "# ", stdout);

    /* read assignment statements in input file */
    Assignment *asmt = Assignment_readEld(ifp);

    Assignment_setInt(asmt, "nbins", &nbins, MANDATORY);
    Assignment_setInt(asmt, "Haploid sample size", &twoNsmp, MANDATORY);
    Assignment_free(asmt);
    asmt = NULL;

    /* fix up interval variables */
    if(outInterval == 0)
        eprintf("Error: output interval (--outInterval) must be > 0");
    if(swapInterval == 0)
        eprintf("Error: output interval (--swapInterval) must be > 0");

    if(outInterval % swapInterval) {
        fprintf(stderr, "Adjusting:"
                " outInterval must be a multiple of swapInterval\n");
        fprintf(stderr, " Old values: outInterval=%u swapInterval=%u\n",
                outInterval, swapInterval);
        outInterval = (1 + outInterval / swapInterval) * swapInterval;
        fprintf(stderr, " New values: outInterval=%u swapInterval=%u\n",
                outInterval, swapInterval);
    }

    if(nreps % swapInterval) {
        fprintf(stderr, "Adjusting:"
                " nreps must be a multiple of swapInterval\n");
        fprintf(stderr, " Old values: nreps=%u swapInterval=%u\n",
                nreps, swapInterval);
        nreps = (1 + nreps / swapInterval) * swapInterval;
        fprintf(stderr, " New values: nreps=%u swapInterval=%u\n",
                nreps, swapInterval);
    }

    printf("# %-35s = %lg\n", "mutation rate per nucleotide", u);
    printf("# %-35s = %d\n", "number of parallel chains", nchains);
    printf("# %-35s = %lg\n", "proposal width: 1/2N", dNinv);
    printf("# %-35s = %lg\n", "proposal width: t", dt);
    printf("# %-35s = %lg\n", "temperature", temp);
    printf("# %-35s = %d\n", "Haploid sample size", twoNsmp);
    printf("# %-35s = %u\n", "MCMC Iterations", nreps);
    printf("# %-35s = %u\n", "output interval", outInterval);
    printf("# %-35s = %u\n", "swap interval", swapInterval);

    cc = (double *) malloc(nbins * sizeof(cc[0]));
    checkmem(cc, __FILE__, __LINE__);

    sigdsq_obs = (double *) malloc(nbins * sizeof(sigdsq_obs[0]));
    checkmem(sigdsq_obs, __FILE__, __LINE__);

    rval = read_data(ifp, nbins, twoNsmp, cc, sigdsq_obs);
    if(rval != nbins)
        eprintf("ERR@%s:%d: Couldn't read %d lines of data from \"%s\"."
                " Only found %d lines", nbins, ifname, rval);

    /* convert centimorgans to recombination rates */
    for(i = 0; i < nbins; ++i)
        cc[i] *= 0.01;

    printf("# %-35s = %s\n", "Model", method);

    Model      *model = Model_alloc(method, twoNsmp);

    checkmem(model, __FILE__, __LINE__);

    pthread_t   threads[nchains];
    Chain      *ch, *chain;

    chain = Chain_new(0, nchains, nreps, u, dt, dNinv, temp, nbins,
                      sigdsq_obs, cc, ph_init, model, odeAbsTol,
                      odeRelTol, swapInterval);

    Chain_sanityCheck(chain, __FILE__, __LINE__);
    /*Chain_prStateAddr(chain, 0, __FILE__,__LINE__); */
    fflush(stdout);
    fprintf(stderr, "Launching %d threads...\n", nchains);

    i = 0;
    for(ch = chain; ch != NULL; ch = Chain_next(ch)) {
        rval = pthread_create(&threads[i], NULL, runChain, ch);
        if(rval)
            eprintf("ERR@%s:%d: pthread_create returned %d\n",
                    __FILE__, __LINE__, rval);
        ++i;
    }

    fflush(stdout);
    fprintf(stderr, "%s: Beginning %u reps...\n", __func__, nreps);

    Chain_sanityCheck(chain, __FILE__, __LINE__);
    Chain_printHdr(chain, stdout);
    for(rep = swapInterval; rep < nreps; rep += swapInterval) {

        PRSTAT("waiting for data");
        Chain_waitForData(chain);
        PRSTAT("got data");
        if(rep % outInterval == 0)
            Chain_printState(chain, stdout);
        Chain_setDataNeeded(chain);
        Chain_unlock(chain);
        Chain_signal(chain);
        fflush(stdout);
    }
    printf("#Finished %u iterations\n", nreps);
    Chain_sanityCheck(chain, __FILE__, __LINE__);

    /* wait for threads to finish */
    for(i = 0; i < nchains; ++i) {
        rval = pthread_join(threads[i], NULL);
        if(rval)
            eprintf("ERR@%s:%d: pthread_join returned %d\n",
                    __FILE__, __LINE__, rval);
        fprintf(stderr, " %2d threads have finished\n", i + 1);
    }

    fprintf(stderr, "Back from threads\n");

    PopHist    *bestPh = PopHist_newEmpty(PopHist_nepoch(ph_init));
    double      bestLnObj;

    Chain_bestFit(chain, &bestLnObj, bestPh);
    Chain_sanityCheck(chain, __FILE__, __LINE__);
    printf("# Optimum: lnObj=%g\n", bestLnObj);
    PopHist_print_comment(bestPh, "#    ", stdout);

    for(ch = chain; ch != NULL; ch = Chain_next(ch)) {
        unsigned    naccpt = Chain_naccepted(ch);
        unsigned    nswap = Chain_nswapped(ch);

        fprintf(stderr, "Chain %u: accepted %u/%u=%lf;"
                " swapped %u/%u=%lf\n",
                Chain_which(ch),
                naccpt,
                nreps,
                naccpt / (double) nreps,
                nswap, nreps, nswap / (double) nreps);
        Chain_sanityCheck(chain, __FILE__, __LINE__);
    }

    Chain_sanityCheck(chain, __FILE__, __LINE__);
    Chain_free(chain);
    if(ifname)
        free(ifname);
    EpochLink_free(linkedList);
    PopHist_free(ph_init);
    PopHist_free(bestPh);
    Model_free(model);
    free(cc);
    free(sigdsq_obs);

    pthread_exit(NULL);
    return 0;
}
