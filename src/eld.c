/**
@file eld.c
@page eld
@brief Estimate LD from genetic data.

`eld`, a program that estimates LD from genetic data
==================================================================

Parameter values, including those describing population history, may
be read either from the initialization file `ldpsiz.ini` or specified
on the command line. Command-line arguments override values in the
initialization file. 

`eld` reads data in "gtp" format and estimates two measures of LD:
\f$\hat\sigma_d^2\f$ and \f$r^2\f$. The program compares all pairs of
sites within a window that slides across the data set. These pairs are
tabulated into bins based on the distance (in cM) that separates
them. For each bin, the program estimates and prints the two measures
of LD.  The maximal separation between pairs of SNPs is determined by
the `--window` (or `-w`) option.

Optionally, the program calculates a moving-blocks bootstrap for
\f$\hat\sigma_d^2\f$. To turn on this option, use the `--bootreps`
option discussed below. The program will then print lower and upper
confidence limits along with each estimate of \f$\hat\sigma_d^2\f$.
By default, the bootstrap includes 95% of the sampling distribution
distribution. This fraction can be modified using the initialization
file. Optionally, the bootstrap may be stored in a separate file for
use by other programs. To accomplish this, use the command-line option 
`--bootfile` or set `bootfile` in `ldpsiz.ini`.

By default, the program uses as many threads as there are cores on the
machine. To change this, use the command-line option `--threads` or
set `nthreads` in the file `ldpsiz.ini`.

    usage: eld [options] input_file_name
       where options may include:
       -b \<x\> or --blocksize \<x\>
          SNPs per bootstrap block
       -h or --help
          print this message
       -i \<x\> or --interval \<x\>
          for debugging: increase speed by examining every x'th focal SNP
       -n \<x\> or --nbins \<x\>
          tablulate values into x bins
       -r \<x\> or --bootreps \<x\>
          number of bootstrap replicates (def: 0)
       -t \<x\> or --threads \<x\>
          number of threads (default is auto)
       -v     or --verbose
          more output
       -w \<x\> or --window \<x\>
          window size in cM

@copyright Copyright (c) 2014, Alan R. Rogers 
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "array.h"
#include "assign.h"
#include "boot.h"
#include "fileindex.h"
#include "ini.h"
#include "misc.h"
#include "readgtp.h"
#include "spectab.h"
#include "tabulation.h"
#include "window.h"
#include <assert.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/** Data required by a single thread */
typedef struct ThreadArg {
    unsigned    thisThread;              /**< index of thread */
    int         nbins;                   /**< number of tabulation bins*/
    int         ploidy;                  /**< 1=haploid; 2=diploid */
    int         overflow;                /**< 0=no overflow */
    unsigned    nGtype;                  /**< number of genotypes in data */

    /**
     * For production, set sampling_interval=1.
     * For debugging, set to (say) 5 to process only every 5th
     * SNP. This makes execution faster. 
     */
    long        sampling_interval;

    /**
     * Maximum distance in centimorgans between pairs of SNPs used in
     * LD calculations.
     */
    double      window_cm;

    const char *ifname;        /**< input file name */

    /**
     * In a multithreaded environment, each thread operates only on a
     * portion of the chromosome. The Threadbounds object, tb,
     * specifies the SNPs to be considered by the current thread.
     */
    ThreadBounds *tb;

    /**
     * Holds tabulated values of numerator and denominator of sigdsq
     * and also of rsq.
     */
    Tabulation *tab;

    /// Holds data for site frequency spectrum.
    Spectab *spectab;

    /**
     * Data used in bootstrap. Not locally owned.
     */
    Boot       *boot;

    /**
     * An index into the input file, which allows the current thread
     * to find the SNPs it needs to analyze.
     */
    FileIndex  *fndx;
} ThreadArg;

static void usage(void);
void        ThreadArg_print(ThreadArg * targ, FILE * ofp);
void       *threadfun(void *varg);

/*pthread_mutex_t mutex_stdout;*/

/** Print ThreadArg object */
void ThreadArg_print(ThreadArg * targ, FILE * ofp) {
    fprintf(ofp, "ThreadArg %d\n", targ->thisThread);
    fprintf(ofp, "  ifname=\"%s\"\n", targ->ifname);
    fprintf(ofp, "  ploidy=%d\n", targ->ploidy);
    fprintf(ofp, "  sampling_interval=%ld\n", targ->sampling_interval);
    fprintf(ofp, "  nbins=%d\n", targ->nbins);
    fprintf(ofp, "  window_cm=%lg\n", targ->window_cm);
    ThreadBounds_print(targ->tb, 1, ofp);
    Tabulation_print(targ->tab, ofp);
    Spectab_print(targ->spectab, ofp);
}

/**
 * `eld` runs a copy of this function within each thread. The
 * function opens the input file, which gives the thread its own input
 * buffer so that it can move around in the file without locking it.
 * Then it examines pairs of SNPs within a window that slides across
 * the SNPs it is responsible for, placing the results into an object
 * of type Tabulation.
 *
 * @param varg A void pointer to an object of type ThreadArg.
 */
void       *threadfun(void *varg) {
    long        i;
    ThreadArg  *targ = (ThreadArg *) varg;

    FILE       *ifp = fopen(targ->ifname, "r");

    if(ifp == NULL) {
        fprintf(stderr, "threadfun: can't open file \"%s\"\n", targ->ifname);
        pthread_exit(NULL);
    }

    /* set up sliding window */
    Window     *window = Window_new(targ->window_cm,
                                    ifp,
                                    targ->sampling_interval,
                                    targ->ploidy);

    /* populate window */
    fseek(ifp, ThreadBounds_seekpos_initial(targ->tb), SEEK_SET);
    for(i = ThreadBounds_ndx_initial(targ->tb);
        i < ThreadBounds_ndx_1stFocal(targ->tb); ++i) {
        Window_nextSNP(window, targ->boot);
    }

    /* slide window */
    for(i = ThreadBounds_ndx_1stFocal(targ->tb);
        i <= ThreadBounds_ndx_lastFocal(targ->tb); ++i) {

        int         status = Window_advance(window, targ->tab,
                                            targ->spectab,
                                            targ->boot, i);

        if(Tabulation_overflow(targ->tab)) {
            targ->overflow = 1;
            break;
        }

        if(status == EOF)
            break;

#if 0
        if(i % 10000 == 0)
            fprintf(stderr, "threadfun SNP #%8ld\n", i);
#endif
    }

    /* Give calling function access to the number of genotypes. */
    targ->nGtype = Window_nGtype(window);

    fclose(ifp);
    Window_free(window);
    pthread_exit(NULL);
}

/** Print usage message and abort */
static void usage(void) {
    fprintf(stderr, "usage: eld [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-b <x> or --blocksize <x>", "SNPs per bootstrap block");
    tellopt("-h or --help", "print this message");
    tellopt("-i <x> or --interval <x>",
            "for debugging: increase speed by examining every x'th focal SNP");
    tellopt("-n <x> or --nbins <x>", "tablulate values into x bins");
    tellopt("-R <x> or --recombination <x>",
            "set recombination rate/generation for adjacent sites");
    tellopt("-r <x> or --bootreps <x>",
            "number of bootstrap replicates (def: 0)");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-v     or --verbose", "more output");
    tellopt("-w <x> or --window <x>", "window size in cM");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

    int         optndx;

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"blocksize", required_argument, 0, 'b'},
        {"help", no_argument, 0, 'h'},
        {"interval", required_argument, 0, 'i'},
        {"nbins", required_argument, 0, 'n'},
        {"bootreps", required_argument, 0, 'r'},
        {"threads", required_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"window", required_argument, 0, 'w'},
        {NULL, 0, NULL, 0}
    };

    double      windowsize_cm = 0.3;
    double      confidence = 0.95;
    int         nthreads = 0;   /* number of threads to launch */
    long        sampling_interval = 1;
    long        bootreps = 0;
    long        blocksize = 300;
    int         nbins = 25;
    int         verbose = 0;
    int         ploidy = 1;     /* default is haploid. */
    const int   folded = true;
    unsigned    nGtype, twoNsmp;

    FILE       *ifp = NULL;
    time_t      currtime = time(NULL);
    gsl_rng    *rng;
    char       *ifname = NULL;
    int         i, tndx;
    int         chromosome = -99;
    long        nSNPs = 0;
    Tabulation **tab;
    Spectab   **spectab;
    Boot      **boot = NULL;
    BootConf   *bc = NULL;
    char        bootfilename[FILENAMESIZE] = { '\0' };
    char        simcmd[1000] = { '\0' };
    unsigned    jobid;
    {
        char s[200];
        snprintf(s, sizeof(s), "%u %s", getpid(), ctime(&currtime));
        jobid = hash(s);
    }

    printf("#################################\n");
    printf("# eld: estimate LD and spectrum #\n");
    printf("#################################\n");

    putchar('\n');
#ifdef __TIMESTAMP__
    printf("# Program was compiled: %s\n", __TIMESTAMP__);
#else
    printf("# Program was compiled: %s, %s\n", __TIME__, __DATE__);
#endif
    printf("# Program was run     : %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    /* import definitions from initialization file */
    Ini        *ini = Ini_new(INIFILE);

    if(ini) {
        Ini_setLong(ini, "blocksize", &blocksize, !MANDATORY);
        Ini_setLong(ini, "samplingInterval", &sampling_interval, !MANDATORY);
        Ini_setInt(ini, "nbins", &nbins, !MANDATORY);
        Ini_setDbl(ini, "confidence", &confidence, !MANDATORY);
        Ini_setInt(ini, "nthreads", &nthreads, !MANDATORY);
        Ini_setDbl(ini, "windowCm", &windowsize_cm, !MANDATORY);
        Ini_free(ini);
        ini = NULL;
    }

    /* command line arguments */
    for(;;) {
        i = getopt_long(argc, argv, "b:f:hi:n:R:r:t:vw:W:", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'b':
            blocksize = strtod(optarg, NULL);
			if(blocksize <= 0) {
				fprintf(stderr,
						"%s:%d: bad argument to -b or --blocksize: \"%s\"\n",
						__FILE__,__LINE__,optarg);
				usage();
			}
            break;
        case 'h':
            usage();
            break;
        case 'i':
            sampling_interval = strtol(optarg, NULL, 10);
            break;
        case 'n':
            nbins = strtod(optarg, NULL);
            break;
        case 'r':
            bootreps = strtol(optarg, NULL, 10);
            break;
        case 't':
            nthreads = strtol(optarg, NULL, 10);
            break;
        case 'v':
            verbose = 1;
            break;
        case 'w':
            windowsize_cm = strtod(optarg, NULL);
            break;
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
            eprintf("ERR@%s:%d: couldn't open file \"%s\"\n",
                    __FILE__, __LINE__, ifname);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    myassert(ifname != NULL);
    myassert(ifp != NULL);

    if(bootreps > 0) {
        // Generate bootstrap output file name from input file name.
        // Strip suffix; add -jobid.boot
        char *basename = strrchr(ifname, '/'); // linux only!!!
        if(basename == NULL)
            basename = ifname;
        snprintf(bootfilename, sizeof(bootfilename), "%s", basename);
        char suffix[30];
        snprintf(suffix, sizeof(suffix), "-%x.boot", jobid);
        replaceSuffix(bootfilename, sizeof(bootfilename), suffix,
                      strlen(suffix));
    }

    /* Read assignments from header of gtp file */
    Assignment *asmt = Gtp_readHdr(ifp);

    Assignment_setInt(asmt, "chromosome", &chromosome, !MANDATORY);
    Assignment_setString(asmt, "sim cmd", simcmd, sizeof(simcmd), !MANDATORY);
    Assignment_setInt(asmt, "ploidy", &ploidy, !MANDATORY);
    Assignment_free(asmt);
    asmt = NULL;

    // 1st pass through file: count SNPs and figure out which will be
    // read by each thread.
    fprintf(stderr, "Indexing file \"%s\"...\n", ifname);
    FileIndex  *fndx = FileIndex_readFile(ifp);

    // Number of SNPs and number genotypes per SNP
    nSNPs = FileIndex_nSNPs(fndx);
    nGtype = FileIndex_nGtype(fndx);
    if(nSNPs == 0 || nGtype == 0) {
        FileIndex_printSummary(fndx, stderr);
        eprintf("ERR@%s:%d: Input file has no data\n", __FILE__, __LINE__);
    }
    assert(ploidy == FileIndex_ploidy(fndx));
    twoNsmp = nGtype * ploidy; // haploid sample size

    /* Number of threads */
    if(nthreads == 0)
        nthreads = lround(getNumCores());

    double      min_cm = FileIndex_getMapPos(fndx, 0);
    double      max_cm = FileIndex_getMapPos(fndx,
                                             FileIndex_nSNPs(fndx) - 1);
    double      maxRange = max_cm - min_cm;

    /* If necessary, reduce window to size of chromosome */
    if(windowsize_cm > maxRange) {
        windowsize_cm = maxRange;

        printf("# Warning: Reducing windowsize to sequence length\n");
        printf("#   windowsize_cm: %lg\n", windowsize_cm);
    }

    ThreadBounds *tb = ThreadBounds_new(nthreads,
                                        windowsize_cm,
                                        fndx);

    FileIndex_free(fndx);
    fndx = NULL;
    fclose(ifp);

    if(verbose)
        ThreadBounds_print(tb, nthreads, stdout);

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned) currtime);

    // allocate arrays
    DblArray *sigdsq = DblArray_new((unsigned long) nbins);
    checkmem(sigdsq, __FILE__, __LINE__);

    DblArray *rsq = DblArray_new((unsigned long) nbins);
    checkmem(rsq, __FILE__, __LINE__);


    DblArray *separation = DblArray_new((unsigned long) nbins);
    checkmem(separation, __FILE__, __LINE__);

    ULIntArray *nobs = ULIntArray_new((unsigned long) nbins);
    checkmem(nobs, __FILE__, __LINE__);
    
    unsigned spdim = specdim(twoNsmp, folded); 
    ULIntArray *spectrum = ULIntArray_new((unsigned long) spdim);
    checkmem(spectrum, __FILE__, __LINE__);

    tab = (Tabulation **) malloc(nthreads * sizeof(tab[0]));
    checkmem(tab, __FILE__, __LINE__);

    spectab = (Spectab **) malloc(nthreads * sizeof(spectab[0]));
    checkmem(spectab, __FILE__, __LINE__);

    for(i = 0; i < nthreads; ++i) {
        tab[i] = Tabulation_new(windowsize_cm, nbins);
        spectab[i] = Spectab_new(twoNsmp, folded);
    }

    // Initially, all boot pointers are set to NULL. If they remain
    // NULL, then each thread will know not do a
    // bootstrap. Bootstrap is done only if boot[i] is not NULL.
    boot = malloc(nthreads * sizeof(*boot));
    checkmem(boot, __FILE__, __LINE__);
    memset(boot, 0, nthreads * sizeof(*boot));

    myassert(boot != NULL);

    for(i = 0; i < nthreads; ++i)
        myassert(boot[i] == NULL);

    if(bootreps > 0) {

		printf("%s:%d: blocksize=%ld\n",
			   __FILE__,__LINE__,blocksize);fflush(stdout);
        boot[0] = Boot_new(nSNPs, bootreps, twoNsmp, folded, blocksize,
                           windowsize_cm,
                           nbins, rng);
        myassert(boot[0] != NULL);

        for(i = 1; i < nthreads; ++i)
            boot[i] = Boot_dup(boot[0]);
    }

    // targ[i] = ptr to i'th ThreadArg object.
    ThreadArg  *targ = malloc(nthreads * sizeof(targ[0]));

    checkmem(targ, __FILE__, __LINE__);
    for(tndx = 0; tndx < nthreads; ++tndx) {
        targ[tndx].thisThread = tndx;
        targ[tndx].ifname = ifname;
        targ[tndx].ploidy = ploidy;
        targ[tndx].nGtype = 0;
        targ[tndx].sampling_interval = sampling_interval;
        targ[tndx].nbins = nbins;
        targ[tndx].window_cm = windowsize_cm;
        targ[tndx].tab = tab[tndx];
        targ[tndx].spectab = spectab[tndx];
        targ[tndx].boot = boot[tndx];
        targ[tndx].tb = &tb[tndx];
        targ[tndx].overflow = 0;
    }

    // echo parameters
    if(strlen(simcmd) > 0)
        printf("# %s = %s\n", "sim cmd", simcmd);
    if(chromosome >= 0)
        printf("# %-35s = %d\n", "chromosome", chromosome);
    printf("# %-35s = %x\n", "JobId", jobid);
    printf("# %-35s = %lg\n", "Window size (cM)", windowsize_cm);
    printf("# %-35s = %s\n", "Input file", ifname);
    printf("# %-35s = %ld\n", "nSNPs", nSNPs);
    printf("# %-35s = %d\n", "nbins", nbins);
    printf("# %-35s = %ld\n", "sampling interval", sampling_interval);
    printf("# %-35s = %d\n", "nthreads", nthreads);
    printf("# %-35s = %ld\n", "bootstrap replicates", bootreps);
    if(bootfilename[0])
        printf("# %-35s = %s\n", "bootstrap output file", bootfilename);
    else
        printf("# %-35s = %s\n", "bootstrap output file", "none");
    printf("# %-35s = %u\n", "Number of genotypes", nGtype);
    printf("# %-35s = %d\n", "Ploidy", ploidy);
    printf("# %-35s = %d\n", "Haploid sample size", nGtype * ploidy);

    fflush(stdout);

    pthread_t  *thread;

    thread = malloc(nthreads * sizeof(pthread_t));
    checkmem(thread, __FILE__, __LINE__);

    fflush(stdout);
    fprintf(stderr, "Launching %d threads...\n", nthreads);
    for(tndx = 0; tndx < nthreads; ++tndx) {
        i = pthread_create(&thread[tndx], NULL, threadfun, &targ[tndx]);
        if(i)
            eprintf("ERR@%s:%d: pthread_create returned %d\n",
                    __FILE__, __LINE__, i);
    }

    fflush(stdout);
    /* wait for threads to finish */
    for(tndx = 0; tndx < nthreads; ++tndx) {
        void       *status;

        i = pthread_join(thread[tndx], &status);
        if(i)
            eprintf("ERR@%s:%d: pthread_join returned %d\n",
                    __FILE__, __LINE__, i);
        fprintf(stderr, " %2d threads have finished\n", tndx + 1);
    }

    fprintf(stderr, "Back from threads\n");

    assert(nGtype == targ[nthreads - 1].nGtype);

    /* aggregate tabulations from threads */
    for(tndx = 1; tndx < nthreads; ++tndx) {
        Tabulation_plus_equals(tab[0], tab[tndx]);
        Spectab_plus_equals(spectab[0], spectab[tndx]);
    }

    if(Tabulation_report(tab[0], separation, nobs, sigdsq, rsq)) {
        fprintf(stderr, "sigdsq data are invalid because"
                " Tabulation overflowed.\n");
        fprintf(stderr, "   Each bin can hold %lu comparisons.\n", ULONG_MAX);
        fprintf(stderr, "   Use fewer SNPs or increase --nbins.\n");
    }else{
        if(bootreps > 0) {
            /* aggregate boot structures from threads */
            for(tndx = 1; tndx < nthreads; ++tndx)
                Boot_plus_equals(boot[0], boot[tndx]);
            bc = BootConf_new(boot[0], confidence);
            myassert(bc);

            if(bootfilename[0]) {
                FILE       *bootfile = fopen(bootfilename, "w");

                Boot_dump(boot[0], bootfile);
                fclose(bootfile);
            }
        }

        putchar('\n');
        printf("# Estimates of sigdsq\n");
        printf("#%12s: mean centimorgans separating pairs of SNPs\n", "cM");
        printf("#%12s: number of observations\n", "nobs");
        if(bootreps > 0)
            BootConf_printHdr(bc, stdout);
        printf("#%10s %11s %10s", "cM", "sigdsq", "nobs");
        if(bootreps > 0)
            printf(" %10s %10s", "loLD", "hiLD");
        printf(" %11s", "rsq");
        putchar('\n');
        for(i = 0; i < nbins; ++i) {
            // "separation" is is distance in cm.
            printf("%11.8lf %11.8lf %10lu",
                   DblArray_get(separation, i),
                   DblArray_get(sigdsq,i),
                   ULIntArray_get(nobs, i));
            if(bootreps > 0)
                printf(" %10.8lf %10.8lf",
                       BootConf_lowBound(bc, i), BootConf_highBound(bc, i));
            printf(" %11.8lf", DblArray_get(rsq,i));
            putchar('\n');
        }

        long unsigned nSpec = Spectab_report(spectab[0], spectrum);
        assert(nSpec == nSNPs);
        
        putchar('\n');
        printf("# %s site frequency spectrum",
               (folded ? "Folded" : "Unfolded"));
        if(bootreps > 0)
            printf(" with %0.1lf%% confidence bounds", 100*confidence);
        putchar('\n');
        printf("# %lu segregating sites\n", nSpec);
        printf("#%10s %11s", "count", "spectrum");
        if(bootreps > 0)
            printf(" %10s %10s", "loSpec", "hiSpec");
        putchar('\n');
        for(i=0; i < spdim; ++i)  {
            printf("%11d %11lu", i+1, ULIntArray_get(spectrum,i));
            if(bootreps > 0)
                printf(" %10.2lf %10.2lf",
                       BootConf_loSpecBound(bc, i),
                       BootConf_hiSpecBound(bc, i));
            putchar('\n');
        }
    }

    if(bootreps > 0) {
        for(i = 0; i < nthreads; ++i)
            Boot_free(boot[i]);
    }
    if(bc)
        BootConf_free(bc);
    free(boot);
    DblArray_free(sigdsq);
    DblArray_free(rsq);
    DblArray_free(separation);
    ULIntArray_free(nobs);
    ULIntArray_free(spectrum);
    for(i = 0; i < nthreads; ++i) {
        Tabulation_free(tab[i]);
    }
    free(thread);
    free(tab);
    free(targ);
    ThreadBounds_free(tb);
    free(ifname);
    gsl_rng_free(rng);

#ifdef __APPLE__
    putchar(' ');               /* prevents hang on exit under osx */
#endif
    return 0;
}
