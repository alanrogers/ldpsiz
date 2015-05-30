/**
   @file polya.h
   @brief Header for polya.c
   @author Alan R. Rogers
   @internal
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#ifndef MATCOAL_POLYA_H
#define MATCOAL_POLYA_H

#include "typedefs.h"
#include <stdio.h>

Polya *Polya_new(int n);
void Polya_free(Polya *polya);
double Polya_prob(const Polya *polya, int i, int k);
void Polya_print(const Polya *polya, FILE *ofp);

#endif
