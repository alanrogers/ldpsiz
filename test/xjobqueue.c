/**
 * @file xjobqueue.c
 * @author Alan R. Rogers
 * @brief Test jobqueue.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "jobqueue.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

typedef struct {
    double      arg, result;
} TstParam;

int         jobfunc(void *p);

int jobfunc(void *p) {
    TstParam   *param = (TstParam *) p;

    param->result = (param->arg) * (param->arg);

    return 0;
}

int main(int argc, char **argv) {

    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xjobqueue [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xjobqueue [-v]\n");
    }

    int         i, njobs = 5, nthreads = 5;
    TstParam    jobs[njobs];
    JobQueue   *jq = JobQueue_new(nthreads);

    for(i = 0; i < njobs; ++i) {
        jobs[i].arg = i + 1.0;
        jobs[i].result = -99.0;
        JobQueue_addJob(jq, jobfunc, jobs + i);
    }

    JobQueue_waitOnJobs(jq);
    JobQueue_free(jq);

    for(i = 0; i < njobs; ++i) {
        if(verbose) {
            printf("%d: %lg --> %lg\n", i, jobs[i].arg, jobs[i].result);
            fflush(stdout);
        }
        assert(jobs[i].result == (i + 1.0) * (i + 1.0));
    }

    unitTstResult("JobQueue", "OK");
    return 0;
}
