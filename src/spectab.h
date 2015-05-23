/**
 * @file spectab.h
 * @author Alan R. Rogers
 * @brief Header for spectab.c
 * @copyright Copyright (c) 2015, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SPECTAB_H
#  define LDPSIZ_SPECTAB_H

#  include "typedefs.h"
#  include "misc.h"

/**
 * spectab.h. This file defines the Spectab object, which tabulates
 * data for the site frequency spectrum. 
 */
Spectab    *Spectab_new(unsigned nSamp, int folded);
Spectab    *Spectab_dup(Spectab * old);
void        Spectab_free(Spectab * tab);
void        Spectab_plus_equals(Spectab * x, const Spectab * y);
void        Spectab_print(Spectab * tab, FILE * ofp);
long unsigned Spectab_report(const Spectab * tab, unsigned dim, int *folded,
                             long unsigned *spec);
int         Spectab_equals(const Spectab * x, const Spectab * y);
void        Spectab_dump(const Spectab * tab, FILE * ofp);
Spectab    *Spectab_restore(FILE * ifp);
void        Spectab_sanityCheck(Spectab * tab, const char *file, int line);
void        Spectab_record(Spectab * tab, unsigned alleleCount,
                           unsigned wgt);
long unsigned Spectab_nObs(const Spectab * tab);

#  ifndef NDEBUG
void        Spectab_test(int verbose);
#  endif
#endif
