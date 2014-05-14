/**
 * @file tstpread.c
 * @author Alan R. Rogers
 * @brief Test parallel reads of a single file by multiple threads. 
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include "jobqueue.h"
#include "misc.h"

/* parameter of taskfun */
typedef struct {
    int         thisjob;
    char        fname[30];
} TaskPar;

int         taskfunc(void *p);

int taskfunc(void *p) {
    TaskPar    *param = (TaskPar *) p;
    char       *bp, buff[1000];
    int         lineno = 0;

    FILE       *fp = fopen(param->fname, "r");

    printf("thisjob=%d pid=%lu pthread_self=%lu\n",
           param->thisjob, (long unsigned) getpid(),
           (long unsigned) pthread_self());

    while(1) {

        bp = fgets(buff, sizeof(buff), fp);
        if(bp == NULL)
            break;
        usleep(10);
        printf("%d:%2d:%s", param->thisjob, ++lineno, buff);
    }

    fclose(fp);
    return 0;
}

int main(int argc, char **argv) {

    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: tstpread [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: tstpread [-v]\n");
    }

    int         rval, i, njobs = 3, nthreads = 3;

    TaskPar     jobs[njobs];

    for(i = 0; i < njobs; ++i) {
        jobs[i].thisjob = i + 1;
        snprintf(jobs[i].fname, sizeof(jobs[i].fname), "dummy.gtp");
    }

    JobQueue   *jq = JobQueue_new(nthreads);

    for(i = 0; i < njobs; ++i)
        JobQueue_addJob(jq, taskfunc, jobs + i);

    JobQueue_waitOnJobs(jq);
    JobQueue_free(jq);

    return 0;
}
