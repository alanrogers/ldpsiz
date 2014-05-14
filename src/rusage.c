/**
 * @file rusage.c
 * @author Alan R. Rogers
 * @brief Use getrusage to measure use of system resources.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <time.h>
#include "rusage.h"

void print_rusage(const char *lbl, struct rusage *ru) {
    if(lbl)
        printf("%9s", lbl);
    else
        printf("%9s", "");
    printf(" %7.3lf", ru->ru_utime.tv_sec + 1e-6 * ru->ru_utime.tv_usec);
    printf(" %7.4lf", ru->ru_stime.tv_sec + 1e-6 * ru->ru_stime.tv_usec);
    printf(" %7ld", ru->ru_maxrss);
#if 0
    printf(" %7ld", ru->ru_ixrss);
    printf(" %7ld", ru->ru_idrss);
    printf(" %7ld", ru->ru_isrss);
#endif
    printf(" %7ld", ru->ru_minflt);
    printf(" %7ld", ru->ru_majflt);
    printf(" %7ld", ru->ru_nswap);
#if 0
    printf(" %7ld", ru->ru_inblock);
    printf(" %7ld", ru->ru_oublock);
    printf(" %7ld", ru->ru_msgsnd);
    printf(" %7ld", ru->ru_msgrcv);
    printf(" %7ld", ru->ru_nsignals);
#endif
    printf(" %7ld", ru->ru_nvcsw);
    printf(" %7ld", ru->ru_nivcsw);
    putchar('\n');
}

void print_rusage_header(void) {
    printf("%9s", "label");
    printf(" %7s", "utime");
    printf(" %7s", "stime");
    printf(" %7s", "maxrss");
#if 0
    printf(" %7s", "ixrss");
    printf(" %7s", "idrss");
    printf(" %7s", "isrss");
#endif
    printf(" %7s", "minflt");
    printf(" %7s", "majflt");
    printf(" %7s", "nswap");
#if 0
    printf(" %7s", "inblock");
    printf(" %7s", "oublock");
    printf(" %7s", "msgsnd");
    printf(" %7s", "msgrcv");
    printf(" %7s", "nsignal");
#endif
    printf(" %7s", "nvcsw");
    printf(" %7s", "nivcsw");
    putchar('\n');
}

/**
 * Subtract each element of y from the corresponding element of x.
 */
void rusage_subtractFrom(struct rusage *x, struct rusage *y) {
    x->ru_utime.tv_sec -= y->ru_utime.tv_sec;
    x->ru_utime.tv_usec -= y->ru_utime.tv_usec;
    x->ru_stime.tv_sec -= y->ru_stime.tv_sec;
    x->ru_stime.tv_usec -= y->ru_stime.tv_usec;
    x->ru_maxrss -= y->ru_maxrss;
    x->ru_ixrss -= y->ru_ixrss;
    x->ru_idrss -= y->ru_idrss;
    x->ru_isrss -= y->ru_isrss;
    x->ru_minflt -= y->ru_minflt;
    x->ru_majflt -= y->ru_majflt;
    x->ru_nswap -= y->ru_nswap;
    x->ru_inblock -= y->ru_inblock;
    x->ru_oublock -= y->ru_oublock;
    x->ru_msgsnd -= y->ru_msgsnd;
    x->ru_msgrcv -= y->ru_msgrcv;
    x->ru_nsignals -= y->ru_nsignals;
    x->ru_nvcsw -= y->ru_nvcsw;
    x->ru_nivcsw -= y->ru_nivcsw;
}

/**
 * Add each element of y from the corresponding element of x.
 */
void rusage_addTo(struct rusage *x, struct rusage *y) {
    x->ru_utime.tv_sec += y->ru_utime.tv_sec;
    x->ru_utime.tv_usec += y->ru_utime.tv_usec;
    x->ru_stime.tv_sec += y->ru_stime.tv_sec;
    x->ru_stime.tv_usec += y->ru_stime.tv_usec;
    x->ru_maxrss += y->ru_maxrss;
    x->ru_ixrss += y->ru_ixrss;
    x->ru_idrss += y->ru_idrss;
    x->ru_isrss += y->ru_isrss;
    x->ru_minflt += y->ru_minflt;
    x->ru_majflt += y->ru_majflt;
    x->ru_nswap += y->ru_nswap;
    x->ru_inblock += y->ru_inblock;
    x->ru_oublock += y->ru_oublock;
    x->ru_msgsnd += y->ru_msgsnd;
    x->ru_msgrcv += y->ru_msgrcv;
    x->ru_nsignals += y->ru_nsignals;
    x->ru_nvcsw += y->ru_nvcsw;
    x->ru_nivcsw += y->ru_nivcsw;
}

void tms_print_header(void) {
    printf("%9s", "label");
    printf(" %12s", "utime");
    printf(" %12s", "stime");
    printf(" %12s", "cutime");
    printf(" %12s", "cstime");
    putchar('\n');
}

void tms_print(const char *lbl, struct tms *t) {
    if(lbl)
        printf("%9s", lbl);
    else
        printf("%9s", "");

    /*
     * According to the man page, units are CLOCKS_PER_SEC. On the
     * mac, however, this gives answers that are smaller than those of
     * get_rusage by a factor of 1e4. I am therefore multiplying by
     * this factor.
     */
    double      scale = 1e4 / CLOCKS_PER_SEC;

    printf(" %12.9lf", t->tms_utime * scale);
    printf(" %12.9lf", t->tms_stime * scale);
    printf(" %12.9lf", t->tms_cutime * scale);
    printf(" %12.9lf", t->tms_cstime * scale);
    putchar('\n');
}

/**
 * Subtract each element of y from the corresponding element of x,
 * provided that the result would be >= zero. Otherwise, set the
 * element to zero.
 */
void tms_subtractFrom(struct tms *x, struct tms *y) {
    x->tms_utime = (y->tms_utime > x->tms_utime
                    ? y->tms_utime - x->tms_utime : 0);

    x->tms_stime = (y->tms_stime > x->tms_stime
                    ? y->tms_stime - x->tms_stime : 0);

    x->tms_cutime = (y->tms_cutime > x->tms_cutime
                     ? y->tms_cutime - x->tms_cutime : 0);

    x->tms_cstime = (y->tms_cstime > x->tms_cstime
                     ? y->tms_cstime - x->tms_cstime : 0);
}
