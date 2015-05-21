/**
 * @file window.h
 * @author Alan R. Rogers
 * @brief Header for window.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_WINDOW_H
#define LDPSIZ_WINDOW_H

#include "snp.h"
#include "typedefs.h"
#include "boot.h"
#include "em.h"
#include "misc.h"
#include <stdio.h>

/**
 * Estimate the linkage disequilibrium (LD) between single nucleotide
 * polymorphisms (SNPs) as a function of the map distance between
 * those SNPs. Definition is provided openly here (rather than hidden
 * in window.c) so that functions can be inlined. 
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

int         Window_advance(Window * window, Tabulation * tab, Boot * boot,
                           long lineno);
SNP        *Window_currSNP(Window * window);
void        Window_free(Window * window);
int         Window_nextSNP(Window * window, Boot * boot);
Window     *Window_new(double width_cm, FILE * ifp,
                       long sampling_interval, unsigned ploidy);
unsigned    Window_nGtype(const Window * window);
long        Window_nSNPsRead(const Window * window);

#  ifndef NDEBUG
void        Window_test(int verbose);
#  endif

#endif
