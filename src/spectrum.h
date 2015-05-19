/**
 * @file spectrum.h
 * @author Alan R. Rogers
 * @brief Header for spectrum.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_SPECTRUM_H
#define LDPSIZ_SPECTRUM_H

#include "typedefs.h"

Spectrum *Spectrum_new(unsigned nSamples);
void        Spectrum_free(Spectrum * spec);
void        Spectrum_print(Spectrum * spec, FILE * ofp);
double      Spectrum_sqrDiff(Spectrum *spec1, Spectrum *spec1);
double      Spectrum_value(Spectrum *s, unsigned i);
#endif
