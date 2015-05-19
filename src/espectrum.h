/**
   @file espectrum.h
   @brief Header for espectrum.c.
   @internal
 
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#ifndef MATCOAL_ESPECTRUM_H
#define MATCOAL_ESPECTRUM_H

#include "typedefs.h"

ESpectrum  *ESpectrum_new(unsigned nSamples, PopHist * ph, double errTol);
void        ESpectrum_free(ESpectrum * spectrum);
double      ESpectrum_unfolded(ESpectrum * spectrum, unsigned i);
double      ESpectrum_folded(ESpectrum * spectrum, unsigned i);
unsigned    ESpectrum_nSamples(const ESpectrum * spectrum);

#endif
