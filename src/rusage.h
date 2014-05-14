/**
 * @file rusage.c
 * @author Alan R. Rogers
 * @brief Header for rusage.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef ARR_RUSAGE
#define ARR_RUSAGE

#include <sys/time.h>           /* for rusage */
#include <sys/resource.h>       /* for rusage */
#include <sys/times.h>          /* for tms    */

void        print_rusage_header(void);
void        print_rusage(const char *lbl, struct rusage *ru);
void        rusage_subtractFrom(struct rusage *x, struct rusage *y);
void        tms_print_header(void);
void        tms_print(const char *lbl, struct tms *t);
void        tms_subtractFrom(struct tms *x, struct tms *y);
void        rusage_addTo(struct rusage *x, struct rusage *y);

#endif
