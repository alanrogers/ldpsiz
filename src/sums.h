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

double      dotprod(unsigned n, double x[n], double y[n]);
unsigned    dotprod_int(unsigned n, unsigned x[n], unsigned y[n]);
unsigned    dotprod_char(unsigned n, unsigned char x[n], unsigned char y[n]);
unsigned    sum_and_int(unsigned n, unsigned x[n], unsigned y[n]);
unsigned    sum_and_char(unsigned n, unsigned char x[n], unsigned char y[n]);
unsigned long sum_long(unsigned n, unsigned long x[n]);
unsigned    sum_int(unsigned n, unsigned x[n]);
double      sum_double(unsigned n, double x[n]);
unsigned    sum_char(unsigned n, const unsigned char x[n]);
unsigned    dotprodDiploid(unsigned n, unsigned char x[n], unsigned char y[n]);
unsigned    sumDiploid(unsigned n, const unsigned char x[n]);

void        trickOptimizer(void);
unsigned long sum_long_slow(unsigned n, unsigned long x[n]);
double      dotprod_slow(unsigned n, double x[n], double y[n]);
unsigned    dotprod_int_slow(unsigned n, unsigned x[n], unsigned y[n]);
unsigned    dotprod_char_slow(unsigned n, unsigned char x[n],
                              unsigned char y[n]);
#endif
