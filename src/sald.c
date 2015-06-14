/**
@file sald.c
@page sald
@brief use simplex-simulated annealing to estimate population history

`sald`: use simplex-simulated annealing to estimate population history
======================================================================

To estimate the parameters describing population history, we need to
find values that provide the best fit to data, which include both
linkage disequilibrium (measured by \f$\sigma_d^2\f$) and the site
frequency spectrum. This involves maximization on a complex surface
with lots of local peaks. (To verify this for yourself, see \ref
mcmcld "mcmcld".) To avoid getting stuck on local peaks, `sald` uses
the simplex version of simulated annealing.

Usage
-----

The input data file should be as produced by \ref eld "eld".

Although simulated annealing works pretty well, I often find that
different runs end up on different peaks. Therefore, `sald` is able to
launch multiple simulated annealing jobs, each from a random starting
point. These can run in parallel, on separate threads. The number of
parallel optimizers is set using the `--nOpt` argument described
below. The 0th optimizer always starts from the population history
parameters specified in `ldpsiz.ini` or on the command line. The other
optimizers start from randomly chosen history parameters.

The program also deals with bootstrap data, as provided by \ref eld
"eld". The various boostrap data sets also run in parallel if your
machine has multiple cores.  By default `sald` does not process
bootstrap replicates: use `--bootfile` if you want it to.

By default `sald` uses as many threads as your machine has cores. This
is not a good idea if you are sharing a machine with other users. Set
the number of threads to some smaller number using `--threads`. On
linux or osx, you can use `top` to figure out how many threads are
actually running. If you launch 20 threads but only 10 run at any
given time, your job will run slower. Stop it and launch again with
`--threads 10`.

Simulated annealing works by beginning with a flattened version of the
objective function. In this flattened version, all the peaks are
smaller, so it is easy for the simplex to move from peak to peak. The
extent of flattening is controlled by a parameter called
"temperature". High temperature implies lots of flattening.  The
annealing algorithm runs for awhile at a high temperature, then lowers
the temperature and runs awhile more. The succession of temperatures
and the number of iterations at each temperature is called the
"annealing schedule". You can change the performance of the algorithm
by adjusting this schedule. See the `--initTmptr`, `--nPerTmptr` and
`--tmptrDecay` arguments.

If a bootstrap file is specified using the -f or --bootfile option,
`sald` will also write a file containing the parameter estimates for
each bootstrap replicate. This file has a name like
<input_file_name>-<jobid>.fboot, where <input_file_name> is the base
name of the input file, and <jobid> is a hexadecimal number intended
to uniquely identify the current job. This .fboot file can be used as
input to the program `tabfboot`.

    usage: sald [options] input_file_name
       where options may include:
       -m <method> or --methods <method>
          specify method "Hill", or "Strobeck", or "Hill,Strobeck"
       -n <x> or --twoNsmp <x>
          haploid sample size
       -t <x> or --threads <x>
          number of threads (default is auto)
       -u <x> or --mutation <x>
          mutation rate/generation
       -v <x> or --verbose
          more output
       -f <x> or --bootfile <x>
          read bootstrap file x
       -c <x> or --confidence <x>
          specify confidence level for CIs of parameters
       --twoN <x>
          haploid pop size to x in current epoch
       -T <x> or --time <x>
          length of current epoch (generations)
       -E or --nextepoch
          move to next earlier epoch
       --nOpt <x>
          optimizers per data set
       --initTmptr <x>
          initial temperature
       --nTmptrs <x>
          number of temperatures
       --tmptrDecay <x>
          ratio of successive temperatures
       -i <x> or --nItr <x>
          total number of iterations
       -h or --help
          print this message

@copyright Copyright (c) 2014, 2015, Alan R. Rogers 
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#undef DEBUG

#undef DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

#define DO_LD 1
#define DO_SPEC 1

#include "annealsched.h"
#include "array.h"
#include "assign.h"
#include "boot.h"
#include "espectrum.h"
#include "hill.h"
#include "ini.h"
#include "jobqueue.h"
#include "misc.h"
#include "model.h"
#include "pophist.h"
#include "polya.h"
#include "sasimplex.h"
#include "spectab.h"
#include "tokenizer.h"
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#if 0
extern long tmr_count;
extern clock_t tmr_cputime;
#endif

extern const int setSimplexVersion;
pthread_mutex_t stdoutLock = PTHREAD_MUTEX_INITIALIZER;

/**
 * The total job is divided into tasks, which are placed on a queue.
 * Each thread takes a task from the queue and executes it. This continues
 * until the queue is empty. The TaskArg structure contains the parameters
 * of a single task.
 */
typedef struct TaskArg {
    unsigned    task;
    unsigned    seed;
    int         nbins;
    unsigned    spdim;
    unsigned    twoNsmp;        // haploid sample size
    size_t      ndim;           // number of dimensions in state space
    double      u;              // mutation rate
    double      tolMatCoal;     // tolerance for MatCoal
    double      ftol, xtol;     // controls convergence
    AnnealSched *sched;         // annealing schedule
    double     *stepsize;       // size of initial simplex
    int         nPerTmptr;      // number of iterations per temp
    int         randomStart;    // whether to initialize from random ph
    int         verbose;
    double     *sigdsq_obs;
    double     *c;
    double     *spectrum_obs;
    double     *loBnd;
    double     *hiBnd;
    double     *hiInit;
    const Polya *polya;         // not locally owned
    ODE        *ode;

    unsigned long nIterations;  // count iterations
    int         status;
    double      simplexSize;
    PopHist    *ph;
    double      cost;
} TaskArg;

/** Parameters of cost function--that which is minimized. */
typedef struct CostPar {
    int         nbins;
    unsigned    spdim;
    unsigned    twoNsmp;
    double      u;
    double      tolMatCoal;
    double     *sigdsq;
    double     *c;
    double     *spectrum;
    ODE        *ode;
    PopHist    *ph;
    const Polya *polya;
    unsigned long nIterations;
} CostPar;

// Structure to pass values to TaskArg_new
typedef struct TaskArgArg {
    unsigned seed;
    int nbins;
    unsigned twoNsmp;
    double u;
    double tolMatCoal;
    double ftol;
    double xtol;
    double *stepsize;
    AnnealSched *sched;
    double *loBnd;
    double *hiBnd;
    double *hiInit;
    double odeAbsTol;
    double odeRelTol;
    int nPerTmptr;
    int verbose;
    double *sigdsq_obs;
    double *c;
    ULIntArray *spectrum;
    Model * model;
    PopHist * ph_init;
    const Polya *polya;
} TaskArgArg;

void        usage(void);
int read_data(FILE * ifp, 
              DblArray *cm,
              DblArray *sigdsq,
              ULIntArray *spectrum);
TaskArg    *TaskArg_new(unsigned task, int randomStart, TaskArgArg *taa);
void        TaskArg_free(TaskArg * targ);
int         taskfun(void *varg);
static double costFun(const gsl_vector *x, void *varg);
void        CostPar_print(CostPar * cp);
TaskArg    *TaskArg_best(TaskArg ** tpvec, int n);
void        prHeader(PopHist * ph);

void CostPar_print(CostPar * cp) {
    int         i;

    printf("CostPar: nbins=%d spdim=%u twoNsmp=%u u=%lg"
           " model=%s nIterations=%lu\n",
           cp->nbins, cp->spdim, cp->twoNsmp, cp->u,
           Model_lbl(ODE_model(cp->ode)), cp->nIterations);
    printf("    %15s %15s\n", "c", "sigdsq");
    for(i = 0; i < cp->nbins; ++i)
        printf("    %15.8lg %15.8lg\n", cp->c[i], cp->sigdsq[i]);

    putchar ('\n');
    printf("    %15s %15s\n", "i", "spectrum");
    for(i = 0; i < cp->spdim; ++i)
        printf("    %15d %15.8lg\n", i+1, cp->spectrum[i]);
    PopHist_print_comment(cp->ph, "    ", stdout);
}

TaskArg    *TaskArg_new(unsigned task, int randomStart, TaskArgArg *taa) {
    TaskArg    *targ = malloc(sizeof(TaskArg));

    checkmem(targ, __FILE__, __LINE__);
    assert(targ != NULL);

    targ->ndim = PopHist_nParams(taa->ph_init);

    targ->nIterations = 0;
    targ->nbins = taa->nbins;
    targ->spdim = ULIntArray_dim(taa->spectrum);
    targ->twoNsmp = taa->twoNsmp;
    targ->u = taa->u;
    targ->task = task;
    /* each task gets different seed */
    targ->seed = (taa->seed + (unsigned long long) task) % UINT_MAX;  

    targ->stepsize = memdup(taa->stepsize, targ->ndim * sizeof(targ->stepsize[0]));
    checkmem(targ->stepsize, __FILE__, __LINE__);

    targ->sched = AnnealSched_copy(taa->sched);

    targ->loBnd = memdup(taa->loBnd, targ->ndim * sizeof(targ->loBnd[0]));
    checkmem(targ->loBnd, __FILE__, __LINE__);

    targ->hiBnd = memdup(taa->hiBnd, targ->ndim * sizeof(targ->hiBnd[0]));
    checkmem(targ->hiBnd, __FILE__, __LINE__);

    targ->hiInit = memdup(taa->hiInit, targ->ndim * sizeof(targ->hiInit[0]));
    checkmem(targ->hiInit, __FILE__, __LINE__);

    targ->nPerTmptr = taa->nPerTmptr;
    targ->verbose = taa->verbose;
    targ->cost = -1.0;
    targ->tolMatCoal = taa->tolMatCoal;
    targ->ftol = taa->ftol;
    targ->xtol = taa->xtol;
    targ->ode = ODE_new(taa->model, taa->odeAbsTol, taa->odeRelTol);
    targ->ph = PopHist_dup(taa->ph_init);
    targ->polya = taa->polya; // not duplicated
    targ->randomStart = randomStart;
    targ->status = 0;
    targ->simplexSize = DBL_MAX;

    targ->sigdsq_obs = memdup(taa->sigdsq_obs, taa->nbins * sizeof(targ->sigdsq_obs[0]));
    checkmem(targ->sigdsq_obs, __FILE__, __LINE__);

    // Normalize spectrum as array of doubles
    unsigned i, spdim = ULIntArray_dim(taa->spectrum);
    double spec[spdim];
    unsigned long spsum=0.0;
    for(i=0; i<spdim; ++i)
        spsum += ULIntArray_get(taa->spectrum, i);
    for(i=0; i<spdim; ++i)
        spec[i] = ULIntArray_get(taa->spectrum, i) / ((double) spsum);

    targ->spectrum_obs = memdup(spec, spdim * sizeof(spec[0]));
    checkmem(targ->spectrum_obs, __FILE__, __LINE__);

    targ->c = memdup(taa->c, taa->nbins * sizeof(targ->c[0]));
    checkmem(targ->c, __FILE__, __LINE__);

    return (targ);
}

void TaskArg_free(TaskArg * targ) {
    // polya isn't freed here because it isn't owned locally.
    PopHist_free(targ->ph);
    free(targ->sigdsq_obs);
    free(targ->c);
    free(targ->stepsize);
    free(targ->loBnd);
    free(targ->hiBnd);
    free(targ->hiInit);
    free(targ->spectrum_obs);
    ODE_free(targ->ode);
    AnnealSched_free(targ->sched);
    free(targ);
}

/**
 * Given a vector of TaskArg pointers, return a pointer to the one
 * with lowest cost among those that have converged.
 */
TaskArg    *TaskArg_best(TaskArg ** tpvec, int n) {
    double      bestCost = DBL_MAX;
    int         i, bestNdx = -1;

    for(i = 0; i < n; ++i) {
        if(tpvec[i]->status==GSL_SUCCESS && tpvec[i]->cost < bestCost) {
            bestNdx = i;
            bestCost = tpvec[i]->cost;
        }
    }
    if(bestNdx == -1)
        return NULL;
    return tpvec[bestNdx];
}

void usage(void) {
    fprintf(stderr, "usage: sald [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");

    /* misc */
    tellopt("-m <method> or --methods <method>",
            "specify method \"Hill\", or \"Strobeck\", or \"Hill,Strobeck\"");
    tellopt("-n <x> or --twoNsmp <x>", "haploid sample size");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-u <x> or --mutation <x>", "mutation rate/generation");
    tellopt("-v or --verbose", "more output");
	tellopt("-s <x> or --truncSpec <x>",
			"truncate <x> entries from site frequency spectrum (def: 0)");

    /* bootstrap */
    tellopt("-f <x> or --bootfile <x>", "read bootstrap file x");
    tellopt("-c <x> or --confidence <x>",
            "specify confidence level for CIs of parameters");

    /* population history */
    tellopt("--twoN <x>", "haploid pop size to x in current epoch");
    tellopt("-T <x> or --time <x>", "length of current epoch (generations)");
    tellopt("-E or --nextepoch", "move to next earlier epoch");

    /* number of parallel optimizers */
    tellopt("--nOpt <x>", "optimizers per data set");

    /* annealing */
    tellopt("--initTmptr <x>", "initial temperature");
    tellopt("--nTmptrs <x>", "number of temperatures");
    tellopt("--tmptrDecay <x>", "ratio of successive temperatures");
    tellopt("-i <x> or --nItr <x>", "total number of iterations");

    tellopt("-h or --help", "print this message");
    exit(1);
}

enum inputState {in_header, in_LD, in_spectrum};

/**
 * Read the data file produced by eld.
 *
 * @param[in] ifp Points to file produced by eld.
 * @param[in] nbins The length of arrays involving LD.
 * @param[out] cm An array giving the average separation (in centimorgans)
 * between pairs of SNPs within the various bins.
 * @param[out] sigdsq An array of estimates of sigma_d^2.
 * @param[out] spectrum An array for the site frequency spectrum.
 *
 * @returns number of lines read
 */
int read_data(FILE * ifp, 
              DblArray *cm,
              DblArray *sigdsq,
              ULIntArray *spectrum) {
    char        buff[200];
    int         ntokens, i, j, tokensExpected = 0;
    enum inputState state = in_header;
    Tokenizer  *tkz = Tokenizer_new(50);

    assert(DblArray_dim(cm) == DblArray_dim(sigdsq));
    unsigned long nbins = DblArray_dim(sigdsq);
           

    rewind(ifp);
    /* skip until beginning of data */
    while(state == in_header
          && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        Tokenizer_split(tkz, buff, " \t");  /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens < 2)
            continue;
        if(strcmp(Tokenizer_token(tkz, 0), "cM") == 0) {
            state = in_LD;
            switch (ntokens) {
            case 4: // fall through
            case 6:
                tokensExpected = ntokens;
                break;
            default:
                fprintf(stderr, "Current tokens:");
                Tokenizer_print(tkz, stderr);
                eprintf("ERR@%s:%d: got %d tokens rather than 4 or 6"
                        " (LD header)",
                        __FILE__, __LINE__, ntokens);
            }
        }
    }

    if(state != in_LD)
        die("Couldn't find LD data in input file", __FILE__, __LINE__);

    // read LD data
    i = 0;
    while(state==in_LD
          && i < nbins
          && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        if(strcomment(buff))
            break;

        ntokens = Tokenizer_split(tkz, buff, " \t");    /* tokenize */
        if(ntokens == 0)
            continue;

        if(ntokens != tokensExpected)
            eprintf("ERR@%s:%d: got %d tokens rather than %d"
                    " (reading LD)",
                    __FILE__, __LINE__, ntokens, tokensExpected);

        DblArray_set(cm, i, strtod(Tokenizer_token(tkz, 0), NULL));
        DblArray_set(sigdsq, i, strtod(Tokenizer_token(tkz, 1), NULL));
        ++i;
    }

    // Skip comments separating LD data from spectrum
    while(fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        if(strempty(buff))
            continue;
        if(!strcomment(buff))
            eprintf("%s:%s:%d: no spectrum in input data",
                    __FILE__, __func__, __LINE__);

        Tokenizer_split(tkz, buff, " \t#");  /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens < 2)
            continue;
        if(strcmp(Tokenizer_token(tkz, 1), "spectrum") == 0) {
            state = in_spectrum;
            switch (ntokens) {
            case 2: // fall through
            case 4:
                tokensExpected = ntokens;
                break;
            default:
                fprintf(stderr, "Current tokens:");
                Tokenizer_print(tkz, stderr);
                eprintf("ERR@%s:%d: got %d tokens rather than 2 or 4"
                        " (spectrum header)",
                        __FILE__, __LINE__, ntokens);
            }
            break;
        }
    }

    // read spectrum
    j = 0;
    while(state==in_spectrum
          && j < ULIntArray_dim(spectrum)
          && fgets(buff, (int) sizeof(buff), ifp) != NULL) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        ntokens = Tokenizer_split(tkz, buff, " \t");    /* tokenize */
        if(ntokens == 0)
            continue;

        if(ntokens != tokensExpected)
            eprintf("ERR@%s:%d: got %d tokens rather than %d"
                    " (reading spectrum)",
                    __FILE__, __LINE__, ntokens, tokensExpected);

#ifndef NDEBUG
        unsigned k = strtod(Tokenizer_token(tkz, 0), NULL);
        assert(k == j+1);
#endif
        unsigned long s = strtoul(Tokenizer_token(tkz, 1), NULL, 10);
        ULIntArray_set(spectrum, j, s);
        ++j;
    }

    if(j != ULIntArray_dim(spectrum)) 
        fprintf(stderr,"%s:%s:%d: WARNING! got %d values of site frequency"
                " spectrum. Was expecting %lu.\n",
                __FILE__, __func__, __LINE__,
                j, ULIntArray_dim(spectrum));

    Tokenizer_free(tkz);
    tkz = NULL;
    return i+j;                   /* return number of lines read */
}

int taskfun(void *varg) {
    TaskArg    *targ = (TaskArg *) varg;
    DPRINTF(("%s:%d:%s: enter task %u thread %lu\n",
             __FILE__, __LINE__, __func__, targ->task,
             (unsigned long) pthread_self()));
    const int   rotate = 0;
    gsl_vector *x = gsl_vector_alloc(targ->ndim);
    gsl_vector_view ss = gsl_vector_view_array(targ->stepsize, targ->ndim);
    int         i, status;
    double      size, temperature;

    CostPar     costPar = {
        .nbins = targ->nbins,
        .spdim = targ->spdim,
        .twoNsmp = targ->twoNsmp,
        .tolMatCoal = targ->tolMatCoal,
        .u = targ->u,
        .sigdsq = targ->sigdsq_obs,
        .c = targ->c,
        .spectrum = targ->spectrum_obs,
        .ode = targ->ode,
        .ph = targ->ph,
        .polya = targ->polya,
        .nIterations = 0
    };

    /* Starting point */
    PopHist_to_vector(x, targ->ph);

    /* Initialize method and iterate */
    gsl_multimin_function minex_func = {
        .n = targ->ndim,
        .f = costFun,
        .params = &costPar
    };
    gsl_multimin_fminimizer *minimizer;
    {
        const gsl_multimin_fminimizer_type *fmType
            = gsl_multimin_fminimizer_sasimplex;
        minimizer = gsl_multimin_fminimizer_alloc(fmType, targ->ndim);
        checkmem(minimizer, __FILE__, __LINE__);
    }

    gsl_multimin_fminimizer_set(minimizer, &minex_func, x, &ss.vector);
    sasimplex_random_seed(minimizer, targ->seed);
    if(targ->verbose) {
        printf("Using minimizer %s.\n",
               gsl_multimin_fminimizer_name(minimizer));
        fflush(stdout);
    }

    /* inequality constraints on parameters */
    gsl_vector_view loBnd =
        gsl_vector_view_array(targ->loBnd, targ->ndim);
    gsl_vector_view hiBnd =
        gsl_vector_view_array(targ->hiBnd, targ->ndim);
    sasimplex_set_bounds(minimizer, &loBnd.vector, &hiBnd.vector, &ss.vector);

    /* for random restarts */
    if(targ->randomStart) {
        gsl_vector_view hiInit =
            gsl_vector_view_array(targ->hiInit, targ->ndim);
        sasimplex_randomize_state(minimizer, rotate,
                                  &loBnd.vector, &hiInit.vector,
                                  &ss.vector);
    }

    for(i=0; i < AnnealSched_size(targ->sched); ++i) {
        temperature = AnnealSched_next(targ->sched);
        status = sasimplex_n_iterations(minimizer,
                                        &size,
                                        targ->ftol,
                                        targ->xtol,
                                        (temperature==0.0
                                         ? 2.0*targ->nPerTmptr
                                         : targ->nPerTmptr),
                                        temperature,
                                        targ->verbose);
        if(targ->verbose) {
            printf("%s:tmptr=%6.4lf hsiz=%.5lf vsiz=%.5lf badness=%.5lf x=",
                   __func__,
                   temperature, size, 
                   sasimplex_vertical_scale(minimizer),
                   minimizer->fval);
            pr_gsl_vector(stdout, "%+7.4le", minimizer->x);
            putchar('\n');
            fflush(stdout);
        }
        if(status != GSL_CONTINUE)
            break;
    }
#if 0
    if(!targ->verbose) {
        fflush(stdout);
        fprintf(stderr,"%s:size=%.5lf vscale=%.5lf badness=%.5lf x=",
               __func__,
               size, sasimplex_vertical_scale(minimizer),
               minimizer->fval);
        pr_gsl_vector(stderr, "%+7.4lf", minimizer->x);
    }
    switch (status) {
    case GSL_SUCCESS:
        printf(" converged\n");
        break;
    case GSL_ETOLF:
        printf(" converged on suboptimal point\n");
        break;
    case GSL_ETOLX:
        printf(" flat objective function\n");
        break;
    case GSL_CONTINUE:
        printf(" no convergence\n");
        break;
    default:
        printf(" unknown status: %d\n", status);
    }
    fflush(stdout);
#endif

    DPRINTF(("%s:%d:%u done minimizing\n", __func__, __LINE__, targ->task));

    targ->nIterations += costPar.nIterations;
    costPar.nIterations = 0;
    targ->status = status;
    vector_to_PopHist(targ->ph, minimizer->x);
    targ->cost = minimizer->fval;
    targ->simplexSize = size;

    if(targ->verbose) {
        pthread_mutex_lock(&stdoutLock);
        printf("# %3s", "");
        for(i = 0; i < PopHist_nParams(targ->ph); ++i) {
            if(i % 2 == 0)
                printf(" %14.6lg", PopHist_paramValue(targ->ph, i));
            else
                printf(" %10.3lf", PopHist_paramValue(targ->ph, i));
        }
        printf(" %11.4lg", targ->cost);
        printf(" %11.4lg", targ->simplexSize);
        printf(" %6d\n", targ->status);
        fflush(stdout);
        pthread_mutex_unlock(&stdoutLock);
    }

    gsl_vector_free(x);
    gsl_multimin_fminimizer_free(minimizer);

    DPRINTF(("%s:%d:%s: exit task %u thread %lu\n",
             __FILE__, __LINE__, __func__, targ->task,
             (unsigned long) pthread_self()));
    return 0;
}

/**
 * Sum of differences between observed and expected sigma_d^2 vectors.
 *
 * @param[in] ph Current population history
 * @param[in] u mutation rate per site per generation
 * @param[in] nbins Number of values in vectors obs and c
 * @param[in] sigdsq Vector of nbins values, the observed values of
 * sigdsq.
 * @param[in] c Vector of nbins values, the recombination rates
 * associated with the values in sigdsq.
 */
static double costFun(const gsl_vector *x, void *varg) {
    DPRINTF(("%s %lu entry\n", __func__,
             (unsigned long) pthread_self()));
    CostPar    *arg = (CostPar *) varg;
    PopHist    *ph = arg->ph;
    int         nbins = arg->nbins;
    unsigned    spdim = arg->spdim;
    double      u = arg->u;
    double     *sigdsq = arg->sigdsq;
    double     *c = arg->c;
    double     *spectrum = arg->spectrum;
    ODE        *ode = arg->ode;
    static const double bigval = 1.0 / DBL_EPSILON;

    double      badness, exp_sigdsq[nbins], exp_spectrum[spdim];
    int         i;

    // DEBUG: make a vector xx = abs(x)
    double xx[x->size];
    for(i=0; i < x->size; ++i)
        xx[i] = fabs(gsl_vector_get(x, i));

    C_array_to_PopHist(ph, x->size, xx);  // DEBUG: using xx instead of x

#if 0
    printf("%s:%d:%s: initial PopHist:\n", __FILE__,__LINE__,__func__);
    fflush(stdout);
    PopHist_print(ph, stdout);
    fflush(stdout);
#endif

    double diff;
    badness = 0.0;

#ifdef DO_LD
    // get vector of expected values of sigdsq
    ODE_ldVec(ode, exp_sigdsq, nbins, c, u, ph);

    for(i = 0; i < nbins; ++i) {
        diff = exp_sigdsq[i] - sigdsq[i];
        badness += diff * diff;
    }
#endif

#ifdef DO_SPEC
    {
        // get expected spectrum
        ESpectrum *espec = ESpectrum_new(arg->twoNsmp, ph,
                                         arg->polya, arg->tolMatCoal);
        for(i=0; i < spdim; ++i)
            exp_spectrum[i] = ESpectrum_folded(espec, i+1);

        for(i = 0; i < spdim; ++i) {
            diff = exp_spectrum[i] - spectrum[i];
            badness += diff * diff;
        }
        ESpectrum_free(espec);
    }
#endif

    if(!isfinite(badness)) {
#ifdef DEBUG
        char        buff[20];

        fprintf(stderr, "%s:%s:%d: badness=%lg for",
                __FILE__, __func__, __LINE__, badness);
        for(i = 0; i < PopHist_nParams(ph); ++i) {
            PopHist_paramName(ph, buff, sizeof(buff), i);
            fprintf(stderr, " %s=%lg", buff, PopHist_paramValue(ph, i));
        }
        putc('\n', stderr);
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
        badness = bigval;
    }
#ifndef NDEBUG
    else
        assert(badness < bigval);
#endif
    ++arg->nIterations;

    DPRINTF(("%s %lu exit\n", __func__,
             (long unsigned) pthread_self()));
    return badness;
}

/** Print header */
void prHeader(PopHist * ph) {
    char        buff[30];
    int         i, nPar = PopHist_nParams(ph);

    printf("# Results from all optimizers\n");
    printf("# %3s", "");
    for(i = 0; i < nPar; ++i) {
        PopHist_paramName(ph, buff, sizeof(buff), i);
        if(i % 2 == 0)
            printf(" %14s", buff);
        else
            printf(" %10s", buff);
    }
    printf(" %11s %11s %8s %6s\n", "badness", "simplexSize", "nItr", "convrg");
    fflush(stdout);
}

int main(int argc, char **argv) {
    int         optndx;
    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"bootfile", required_argument, 0, 'f'},
        {"confidence", required_argument, 0, 'c'},
        {"nextepoch", no_argument, 0, 'E'},
        {"tmptrDecay", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {"initTmptr", required_argument, 0, 'I'},
        {"nItr", required_argument, 0, 'i'},
        {"nOpt", required_argument, 0, 'o'},
        {"methods", required_argument, 0, 'm'},
        {"twoN", required_argument, 0, 'N'},
        {"twoNsmp", required_argument, 0, 'n'},
        {"nTmptrs", required_argument, 0, 'p'},
        {"threads", required_argument, 0, 't'},
        {"time", required_argument, 0, 'T'},
        {"mutation", required_argument, 0, 'u'},
        {"verbose", no_argument, 0, 'v'},
		{"truncSpec", required_argument, 0, 's'},
        {NULL, 0, NULL, 0}
    };

    double      u = 1e-8;
    double      ftol = 5e-5;     /* will be multiplied by nbins */
    double      xtol = 1e-6;
    double      tolMatCoal = 1e-3; // controls accuracy of MatCoal
    double      odeAbsTol = 1e-7;
    double      odeRelTol = 1e-3;
    double      confidence = 0.95;
    double      initTmptr = 2.0;  
    double      tmptrDecay = 1.0;
    int         nItr = 1000;     /* total number of iterations */
    int         nPerTmptr;       /* iterations at each temperature */
    int         nTmptrs = 3;     /* number of temperatures */
    const int   folded = true;   // Folded site frequency spectrum
	int         truncSFS = 0;    // Truncate site frequency spectrum
    double      lo2N = 400.0, hi2N = 1e8, lo2Ninit = 1000.0;
    double      loT = 1.0, hiT = 5e3, hiTinit = 2000.0;
    double     *stepsize;            /* controls size of initial simplex */
    double      durationEps = 500.0;
    double      twoNinvEps = 0.01;

    int         nbins;
    int         verbose = 0;
    int         twoNsmp = 0;    /* number of haploid samples */
    int         nthreads = 0;   /* number of threads to launch */
    int         nDataSets = 1;  /* 1 + nBootReps */
    int         nOpt = 10;      /* parallel optimization threads
                                 * per data set */
    int         nTasks = 0;     /* total number of tasks */
    long        nBootReps = 0;  /* bootstrap replicates */
    int         curr_epoch = 0, epochPending = 1, phSetFromFile = 0;
    double      curr_t = strtod("Inf", 0);
    double      curr_twoN = 1000.0;
    double      badness;

    char        method[20] = { 0 };
    EpochLink  *linkedList = NULL;
    char        bootfilename[FILENAMESIZE] = { '\0' };
    char        fname[FILENAMESIZE] = { '\0' };
    FILE       *ifp = NULL;

    int         i, j, rval;
    int         rndx, nparams = 0, pndx;
    PopHist    *ph_init = NULL;
    Model      *model = NULL;
    time_t      currtime = time(NULL);
    unsigned    baseSeed = currtime % UINT_MAX;

    unsigned    jobid;
    {
        char s[200];
        snprintf(s, sizeof(s), "%u %s", getpid(), ctime(&currtime));
        jobid = hash(s);
    }


    printf("############################################\n"
           "# sald: fit population history to LD using #\n"
           "#       simulated annealing                #\n"
           "############################################\n");
    putchar('\n');
#ifdef __TIMESTAMP__
    printf("# Program was compiled: %s\n", __TIMESTAMP__);
#endif
    printf("# Program was run     : %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    fflush(stdout);

    /* import definitions from initialization file */
    Ini        *ini = Ini_new(INIFILE);

    if(ini) {
        Ini_setDbl(ini, "confidence", &confidence, !MANDATORY);
        Ini_setString(ini, "methods", method, sizeof(method), !MANDATORY);
        Ini_setDbl(ini, "mutation", &u, !MANDATORY);
        Ini_setDbl(ini, "initTmptr", &initTmptr, !MANDATORY);
        Ini_setDbl(ini, "tmptrDecay", &tmptrDecay, !MANDATORY);
        Ini_setInt(ini, "nTmptrs", &nTmptrs, !MANDATORY);
        Ini_setInt(ini, "nthreads", &nthreads, !MANDATORY);
        Ini_setInt(ini, "nPerTmptr", &nPerTmptr, !MANDATORY);
        Ini_setInt(ini, "nOpt", &nOpt, !MANDATORY);
        if(Ini_setEpochLink(ini, &linkedList, !MANDATORY))
            phSetFromFile = 1;
    }
    Ini_free(ini);
    ini = NULL;

    /* command line arguments */
    for(;;) {
        i = getopt_long(argc, argv, "c:d:f:BEhi:m:n:N:s:t:T:u:v",
                        myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'c':
            confidence = strtod(optarg, 0);
            break;
        case 'd':
            tmptrDecay = strtod(optarg, 0);
            break;
        case 'E':
            if(!epochPending) {
                fprintf(stderr, "Arg out of place: -E or --nextEpoch\n");
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
            break;
        case 'f':
            snprintf(bootfilename, sizeof(bootfilename), "%s", optarg);
            break;
        case 'I':
            assert(optarg);
            initTmptr = strtod(optarg, 0);
            break;
        case 'm':
            snprintf(method, sizeof(method), "%s", optarg);
            break;
        case 'N':
            myassert(optarg);
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                linkedList = NULL;
                phSetFromFile = 0;
            }
            curr_twoN = strtod(optarg, 0);
            epochPending = 1;
            break;
        case 'n':
            twoNsmp = strtol(optarg, NULL, 10);
            break;
        case 'i':
            nItr = strtol(optarg, NULL, 10);
            break;
        case 'o':
            nOpt = strtol(optarg, NULL, 10);
            break;
        case 'p':
            nTmptrs = strtol(optarg, NULL, 10);
            break;
		case 's':
			truncSFS = strtol(optarg, NULL, 10);
            if(truncSFS <= 0) {
                fprintf(stderr,"Arg following -s or --truncSpec"
						" should be positive. Got %s.\n", optarg);
                usage();
            }
			break;
        case 't':
            nthreads = strtol(optarg, NULL, 10);
            break;
        case 'T':
            myassert(optarg);
            if(phSetFromFile) {
                /* override PopHist from file */
                EpochLink_free(linkedList);
                linkedList = NULL;
                phSetFromFile = 0;
            }
            curr_t = strtod(optarg, 0);
            if(curr_t <= 0) {
                fprintf(stderr,"Arg following -T or --time should be positive."
						" Got %s.\n", optarg);
                usage();
            }
            epochPending = 1;
            break;
        case 'u':
            assert(optarg);
            u = strtod(optarg, 0);
            break;
        case 'v':
            verbose = 1;
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
        snprintf(fname, sizeof(fname), "%s", argv[optind]);
        ifp = fopen(fname, "r");
        if(ifp == NULL)
            eprintf("ERR@%s:%d: Couldn't open %s for input",
                    __FILE__, __LINE__, fname);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(fname[0] != '\0');
    assert(ifp != NULL);

    fflush(stdout);

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
    nparams = PopHist_nParams(ph_init);

    printf("# Initial PopHist:\n");
    PopHist_print_comment(ph_init, "# ", stdout);

    stepsize = malloc(nparams * sizeof(double));
    checkmem(stepsize, __FILE__, __LINE__);
    PopHist_setAllTwoNinv(stepsize, nparams, twoNinvEps);
    PopHist_setAllDuration(stepsize, nparams, durationEps);

    double     *loBnd = malloc(nparams * sizeof(loBnd[0]));
    checkmem(loBnd, __FILE__, __LINE__);

    PopHist_setAllTwoNinv(loBnd, nparams, 1.0/hi2N);
    PopHist_setAllDuration(loBnd, nparams, loT);

    double     *hiBnd = malloc(nparams * sizeof(hiBnd[0]));
    checkmem(hiBnd, __FILE__, __LINE__);

    PopHist_setAllTwoNinv(hiBnd, nparams, 1.0/lo2N);
    PopHist_setAllDuration(hiBnd, nparams, hiT);

    double     *hiInit = malloc(nparams * sizeof(hiInit[0]));
    checkmem(hiInit, __FILE__, __LINE__);

    PopHist_setAllTwoNinv(hiInit, nparams, 1.0/lo2Ninit);
    PopHist_setAllDuration(hiInit, nparams, hiTinit);

    /* read assignment statements in input file */
    Assignment *asmt = Assignment_readEld(ifp);

    Assignment_setInt(asmt, "nbins", &nbins, MANDATORY);
    Assignment_setInt(asmt, "Haploid sample size", &twoNsmp, MANDATORY);
    Assignment_free(asmt);
    asmt = NULL;

    /* scale ftol to number of bins */
    ftol *= nbins;

    /* determine number of iterations at each temperature */
    nPerTmptr = round(nItr/ (double) nTmptrs);
    if(nPerTmptr <= 0)
        nPerTmptr = 1;
    nItr = nPerTmptr * nTmptrs;

    model = Model_alloc(method, twoNsmp);
    AnnealSched *sched = AnnealSched_alloc(nTmptrs, initTmptr, tmptrDecay);

    fflush(stdout);
    // dimension of spectrum
    unsigned spdim = specdim((unsigned) twoNsmp, folded);

    Polya *polya = Polya_new(twoNsmp);

    printf("# %-35s = %x\n", "JobId", jobid);
    printf("# %-35s = %s\n", "Model", method);
    printf("# %-35s = %lg\n", "mutation rate per nucleotide", u);
    printf("# %-35s = %lg\n", "ftol", ftol);
    printf("# %-35s = %lg\n", "xtol", xtol);
    printf("# %-35s = %lg\n", "tolMatCoal", tolMatCoal);
    printf("# %-35s = %lg\n", "odeAbsTol", odeAbsTol);
    printf("# %-35s = %lg\n", "odeRelTol", odeRelTol);
    printf("# %-35s = %s\n", "Fitting LD", (DO_LD ? "true" : "false"));
    printf("# %-35s = %s\n", "Fitting site frequency spectrum",
           (DO_SPEC ? "true" : "false"));

    printf("# %-35s = %lg\n", "Initial temp", initTmptr);
    printf("# %-35s = %lg\n", "final temp", 0.0);
    printf("# %-35s = %lg\n", "tmptr decay", tmptrDecay);
    printf("# %-35s = %d\n", "nTmptrs", nTmptrs);
    printf("# %-35s = %d\n", "iterations in total", nItr);
    printf("# %-35s = %d\n", "iterations/temp", nPerTmptr);
    printf("# %-35s = %d\n", "optimizers/data set", nOpt);
    printf("# %-35s = %d\n", "haploid sample size", twoNsmp);
    printf("# %-35s = [%lg, %lg]\n", "twoN range", lo2N, hi2N);
    printf("# %-35s = [%lg, %lg]\n", "t range", loT, hiT);
    printf("# %-35s = [", "scale of initial simplex");
    for(i = 0; i + 1 < nparams; ++i)
        printf("%lg, ", stepsize[i]);
    printf("%lg]\n", stepsize[i]);

    fflush(stdout);

    DblArray *cc = DblArray_new(nbins);
    checkmem(cc, __FILE__, __LINE__);

    DblArray *sigdsq_obs = DblArray_new(nbins);
    checkmem(sigdsq_obs, __FILE__, __LINE__);

    ULIntArray *spectrum_obs = ULIntArray_new(spdim);
    checkmem(spectrum_obs, __FILE__, __LINE__);

    rval = read_data(ifp, cc, sigdsq_obs, spectrum_obs);
    if(rval != nbins + spdim)
        eprintf("%s:%s:%d: Expected %d lines of data from \"%s\"."
                " Got %d.",
                __FILE__,__func__,__LINE__,
                nbins, fname, rval);

    long unsigned nSNPs = 0;
    for(i=0; i < ULIntArray_dim(spectrum_obs); ++i)
        nSNPs += ULIntArray_get(spectrum_obs, i);
    printf("# %-35s = %lu\n", "nSNPs", nSNPs);

    // convert centimorgans to recombination rates 
    for(i = 0; i < nbins; ++i) {
        double c = DblArray_get(cc, i);
        DblArray_set(cc, i, c * 0.01);
    }

    Boot       *boot = NULL;
    FILE       *bootfile = NULL;

    if(bootfilename[0]) {
        bootfile = fopen(bootfilename, "r");
        if(bootfile == NULL)
            eprintf("ERR@%s:%d: can't open bootfile \"%s\"\n",
                    __FILE__, __LINE__, bootfilename);
        boot = Boot_restore(bootfile);
        if(boot == NULL)
            eprintf("ERR@%s:%d: can't restore from bootfile \"%s\"\n",
                    __FILE__, __LINE__, bootfilename);
        fclose(bootfile);
        bootfile = NULL;
    }

    fflush(stdout);
    if(boot) {
        nBootReps = Boot_nReps(boot);
        Boot_purge(boot);
        if(nBootReps > Boot_nReps(boot)) {
            fprintf(stderr, "Boot_purge removed %ld reps\n",
                    nBootReps - Boot_nReps(boot));
            nBootReps = Boot_nReps(boot);
        }
        assert(nBootReps >= 0);
    } else
        nBootReps = 0;
    printf("# %-35s = %ld\n", "Number of bootstrap replicates", nBootReps);

    fflush(stdout);

    nDataSets = 1 + nBootReps;
    nTasks = nDataSets * nOpt;

    // taskarg[i][j] points to the j'th optimizer on the i'th data
    // set, where the observed data are data set 0. i runs from 0
    // through nDataSets-1, and j from 0 through nOpt-1.
    TaskArg  ***taskarg = malloc(nDataSets * sizeof(taskarg[0]));
    checkmem(taskarg, __FILE__, __LINE__);

    for(i = 0; i < nDataSets; ++i) {
        /* allocate a row of pointers for each data set */
        taskarg[i] = malloc(nOpt * sizeof(taskarg[i][0]));
        checkmem(taskarg[i], __FILE__, __LINE__);
    }

    // These values are all passed to TaskArg_new. Putting them into
    // a structure is less error-prone than a long list of function
    // arguments.  
    TaskArgArg taskargarg = {
        .seed = baseSeed,
        .nbins = nbins,
        .twoNsmp = twoNsmp,
        .u = u,
        .tolMatCoal = tolMatCoal,
        .ftol = ftol,
        .xtol = xtol,
        .stepsize = stepsize,
        .sched = sched,
        .loBnd = loBnd,
        .hiBnd = hiBnd,
        .hiInit = hiInit,
        .odeAbsTol = odeAbsTol,
        .odeRelTol = odeRelTol,
        .nPerTmptr = nPerTmptr,
        .verbose = verbose,
        .sigdsq_obs = DblArray_ptr(sigdsq_obs),
        .c = DblArray_ptr(cc),
        .spectrum = spectrum_obs,
        .model = model,
        .ph_init = ph_init,
        .polya = polya,
    };

    /* create task arguments for the observed sigdsq */
    for(j = 0; j < nOpt; ++j) {
        int randomStart = (j==0 ? true : false);
        taskarg[0][j] = TaskArg_new(j, randomStart, &taskargarg);
    }

    fflush(stdout);

    /* create task arguments for each bootstrap replicate */
    DblArray *sigdsq_curr = DblArray_new(nbins);
    checkmem(sigdsq_curr, __FILE__, __LINE__);
    DblArray *cc_curr = DblArray_new(nbins);
    checkmem(cc_curr, __FILE__, __LINE__);
    ULIntArray *spec_curr = ULIntArray_new(spdim);
    checkmem(spec_curr, __FILE__, __LINE__);

    taskargarg.verbose = false;
    for(rndx = 0; rndx < nBootReps; ++rndx) {

        Boot_get_rep(boot, sigdsq_curr, NULL, cc_curr, NULL,
                     spec_curr, rndx);

        // conv. cM to recombination rate
        for(i = 0; i < nbins; ++i) {
            double cci = DblArray_get(cc_curr, i);
            DblArray_set(cc_curr, i, cci*0.01);
        }

        // These change with each bootstrap replicate
        taskargarg.sigdsq_obs = DblArray_ptr(sigdsq_curr);
        taskargarg.c = DblArray_ptr(cc_curr);
        taskargarg.spectrum = spec_curr;

        for(j = 0; j < nOpt; ++j) {
            int randomStart = (j==0 ? true : false);
            taskarg[rndx + 1][j] = TaskArg_new(j + (1+rndx)*nOpt, randomStart,
                                               &taskargarg);
        }
    }

    fflush(stdout);

    if(nthreads == 0)
        nthreads = getNumCores();

    if(nthreads > nTasks)
        nthreads = nTasks;

    fflush(stdout);
    fprintf(stderr, "Creating %d threads to perform %d tasks\n",
            nthreads, nTasks);

    JobQueue   *jq = JobQueue_new(nthreads);

    if(jq == NULL)
        eprintf("ERR@%s:%d: Bad return from JobQueue_new",
                __FILE__, __LINE__);

    for(i = 0; i < nDataSets; ++i)
        for(j = 0; j < nOpt; ++j)
            JobQueue_addJob(jq, taskfun, taskarg[i][j]);

    if(verbose)
        prHeader(ph_init);
    JobQueue_waitOnJobs(jq);

    fflush(stdout);
    fprintf(stderr, "Back from threads\n");

    TaskArg   **best = malloc(nDataSets * sizeof(best[0]));
    checkmem(best, __FILE__, __LINE__);

    // best[i] points to the best result among all replicate
    // optimizers for data set i.
    for(i = 0; i < nDataSets; ++i)
        best[i] = TaskArg_best(taskarg[i], nOpt);

    if(best[0] == NULL)
        printf("No convergence\n");

    prHeader(ph_init);
    fflush(stdout);
    for(i = 0; i < nOpt; ++i) {
        printf("# %2d:", i);
        for(j = 0; j < PopHist_nParams(taskarg[0][i]->ph); ++j) {
            if(j % 2 == 0)
                printf(" %14.6lg", PopHist_paramValue(taskarg[0][i]->ph, j));
            else
                printf(" %10.3lf", PopHist_paramValue(taskarg[0][i]->ph, j));
        }
        printf(" %11.9lf", taskarg[0][i]->cost);
        printf(" %11.9lf", taskarg[0][i]->simplexSize);
        printf(" %8lu", taskarg[0][i]->nIterations);
        switch(taskarg[0][i]->status) {
        case GSL_SUCCESS:
            printf(" %-9s\n", "Converg");
            break;
        case GSL_ETOLF:
            printf(" %-9s\n", "ETOLF");
            break;
        case GSL_ETOLX:
            printf(" %-9s\n", "ETOLX");
            break;
        case GSL_CONTINUE:
            printf(" %-9s\n","NoConverg");
            break;
        default:
            printf(" %-9d\n", taskarg[0][i]->status);
        }
    }

    fflush(stdout);
    if(best[0] != NULL) {
        badness = best[0]->cost;
        printf("\n# minimum badness = %lg\n", badness);
    }

    char        pname[50];

    if(boot) {
        fprintf(stderr, "Processing bootstrap\n");
        // output with confidence interval
        assert(nBootReps > 0);
        double      low, high;
        char        cibuff[50];
        double     *v = malloc(nBootReps * sizeof(v[0]));

        checkmem(v, __FILE__, __LINE__);

        snprintf(cibuff, sizeof(cibuff), "%2.0lf%%_Conf_Int",
                 round(100 * confidence));
        printf("\n# %5s: %10s %10s %23s\n",
               "param", "initial", "estimated", cibuff);

		fflush(stdout);

        for(pndx = 0; pndx < nparams; ++pndx) {
            long        ngood = 0;
            double      paramVal;

            for(rndx = 0; rndx < nBootReps; ++rndx) {
                if(best[1 + rndx] == NULL)
                    continue;
                paramVal = PopHist_paramValue(best[1 + rndx]->ph, pndx);
                v[ngood++] = paramVal;
            }
            if(ngood == 0) {
                fprintf(stderr,
                        "Warning: All bootstrap reps failed to converge\n");
                break;
            }
            confidenceBounds(&low, &high, confidence, v, ngood);
            rval = PopHist_paramName(ph_init, pname, sizeof(pname), pndx);
            if(rval) {
                fprintf(stderr,
                        "OOPS@%s:%dBad return from PopHist_paramName\n",
                        __FILE__, __LINE__);
            } else if(best[0]) {
                snprintf(cibuff, sizeof(cibuff), "(%5.5lg,%5.5lg)", low,
                         high);
                printf("# %5s: %10.7lg %10.5lg %23s\n", pname,
                       PopHist_paramValue(ph_init, pndx),
                       PopHist_paramValue(best[0]->ph, pndx), cibuff);
            }
			fflush(stdout);
        }
        free(v);
    } else if(best[0]) {
        /* output w/o confidence interval */

        printf("# %10s: %10s %10s\n", "param", "initial", "estimated");
        for(pndx = 0; pndx < nparams; ++pndx) {
            rval = PopHist_paramName(ph_init, pname, sizeof(pname), pndx);
            if(rval)
                fprintf(stderr,
                        "%s:%s:%d:Bad return from PopHist_paramName\n",
                        __func__,__FILE__, __LINE__);
            printf("# %10s: %10.7lg %10.5lg\n",
                   pname,
                   PopHist_paramValue(ph_init, pndx),
                   PopHist_paramValue(best[0]->ph, pndx));
			fflush(stdout);
        }
    }

    fflush(stdout);

    if(best[0]) {
        // find fitted values of sigdsq
        double      sigdsq_fit[nbins];
        ODE        *ode = ODE_new(model, odeAbsTol, odeRelTol);
        ODE_ldVec(ode, sigdsq_fit, nbins,
                  DblArray_ptr(cc), u, best[0]->ph);
        ODE_free(ode);

        printf("\n# Fitted values of sigma_d^2. (cM = centimorgans)\n");
        printf("#%10s %10s\n", "cM", "sigma_d^2");
        for(i = 0; i < nbins; ++i) {
            double cci = DblArray_get(cc, i);
            printf("%11.8lf %10.8lf\n",
                   cci * 100.0,  /*convert to cM*/
                   sigdsq_fit[i]);
        }

        // fitted spectrum
        ESpectrum *espec = ESpectrum_new(twoNsmp, best[0]->ph,
                                         polya, tolMatCoal);
        printf("\n# Fitted site frequency spectrum\n");
        printf("#%5s %10s\n", "i", "spectrum");
        for(i=1; i <= spdim; ++i)
            printf("%6d %10.0lf\n", i,
                   floor(0.5 + nSNPs*ESpectrum_folded(espec, i)));

        ESpectrum_free(espec);
    }

    /* mean and variance of log badness */
    double      m = 0.0, v = 0.0;

    for(i = 0; i < nOpt; ++i) {
        double      lncost = log(taskarg[0][i]->cost);

        m += lncost;
        v += lncost * lncost;
    }
    m /= nOpt;
    v = (v - nOpt * m * m) / (nOpt - 1);
    printf("# log badness:m=%lg sd=%lg\n", m, sqrt(v));

    fflush(stdout);

    // Generate an fboot file, which is a rectangular table. There is
    // one column for each estimated parameter. Row 1 contains the
    // parameter labels. Row 2 is the real data. Each succeeding row
    // contains the estimate from one bootstrap replicate. Replicates
    // that did not converge are omitted.
    if(boot && best[0]) {
        // strip suffix from fname; add suffix -jobid.fboot; open file
        char *basename = strrchr(fname, '/'); // linux only!!!
        if(basename == NULL)
            basename = fname;
        else
            basename += 1;    // skip '/' character
        snprintf(bootfilename, sizeof(bootfilename), "%s", basename);
        char suffix[30];
        snprintf(suffix, sizeof(suffix), "-%x.fboot", jobid);
        replaceSuffix(bootfilename, sizeof(bootfilename), suffix,
                      strlen(suffix));
        bootfile = fopen(bootfilename, "w");
        printf("# %-35s = %s\n", "fboot file name", bootfilename);

        // 1st line of fboot file contains parameter names
        for(pndx = 0; pndx < nparams; ++pndx) {
            rval = PopHist_paramName(ph_init, pname, sizeof(pname), pndx);
            if(rval)
                fprintf(stderr,
                        "%s:%s:%d:Bad return from PopHist_paramName: %d\n",
                        __func__,__FILE__, __LINE__, rval);
            fprintf(bootfile, " %*s", DBL_DIG + 9, pname);
        }
        putc('\n', bootfile); // end of 1st line 

        // each subsequent line contains parameter values
        for(i=0; i <= nBootReps; ++i) {
            if(best[i] == NULL)
                continue;

            // Print parameter values in exponential format using
            // precision DBL_DIG+3.  This is as much as you can do in
            // decimal format.
            for(pndx=0; pndx < nparams; ++pndx)
                fprintf(bootfile, " %.*le", DBL_DIG+3,
                        PopHist_paramValue(best[i]->ph, pndx));
            putc('\n', bootfile);
        }
        fclose(bootfile);
        bootfile = NULL;
    }

    fflush(stdout);

    DblArray_free(sigdsq_curr);
    DblArray_free(cc_curr);
    ULIntArray_free(spec_curr);
    free(best);
    for(i = 0; i < nDataSets; ++i) {
        for(j = 0; j < nOpt; ++j)
            TaskArg_free(taskarg[i][j]);
        free(taskarg[i]);
    }
    free(taskarg);
    free(stepsize);
    EpochLink_free(linkedList);
    PopHist_free(ph_init);
    DblArray_free(cc);
    DblArray_free(sigdsq_obs);
    ULIntArray_free(spectrum_obs);
    if(boot) {
        Boot_free(boot);
    }
    JobQueue_free(jq);
    AnnealSched_free(sched);
    free(loBnd);
    free(hiBnd);
    free(hiInit);
    fprintf(stderr, "sald is finished\n");
	Model_free(model);
    Polya_free(polya);

    return 0;
}
