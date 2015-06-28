/**
   @file tfespectrum.h
   @brief Header for tfespectrum.c.
   @internal
 
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#ifndef MATCOAL_TFESPECTRUM_H
#define MATCOAL_TFESPECTRUM_H

#include "typedefs.h"

TFESpectrum  *TFESpectrum_new(unsigned nSamples, unsigned truncSFS, PopHist *ph,
                         const Polya *polya, double errTol);
void        TFESpectrum_free(TFESpectrum * self);
unsigned    TFESpectrum_dim(const TFESpectrum *self);
double     *TFESpectrum_ptr(TFESpectrum *self);
double      TFESpectrum_atNdx(TFESpectrum * self, unsigned i);
unsigned    TFESpectrum_nSamples(const TFESpectrum * self);
double      TFESpectrum_diff(const TFESpectrum *self, unsigned dim, double s[dim]);
double      TFESpectrum_KLdiverg(const TFESpectrum *self, unsigned dim, double s[dim]);

#endif
