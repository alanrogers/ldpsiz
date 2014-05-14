/**
 * @file sums.h
 * @author Alan R. Rogers
 * @brief Header for sums.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SUMS_H
#define LDPSIZ_SUMS_H

double      dotprod(double *x, double *y, unsigned n);
unsigned    dotprod_int(unsigned *x, unsigned *y, unsigned n);
unsigned    dotprod_char(unsigned char *x, unsigned char *y, unsigned n);
unsigned    sum_and_int(unsigned *x, unsigned *y, unsigned n);
unsigned    sum_and_char(unsigned char *x, unsigned char *y, unsigned n);
unsigned long sum_long(unsigned long *x, unsigned n);
unsigned    sum_int(unsigned *x, unsigned n);
double      sum_double(double *x, unsigned n);
unsigned    sum_char(const unsigned char *x, unsigned n);
unsigned    dotprodDiploid(unsigned char *x, unsigned char *y, unsigned n);
unsigned    sumDiploid(const unsigned char *x, unsigned n);

void        trickOptimizer(void);
unsigned long sum_long_slow(unsigned long *x, unsigned n);
double      dotprod_slow(double *x, double *y, unsigned n);
unsigned    dotprod_int_slow(unsigned *x, unsigned *y, unsigned n);
unsigned    dotprod_char_slow(unsigned char *x, unsigned char *y, unsigned n);

#endif
