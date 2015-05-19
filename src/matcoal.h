/**
   @file matcoal.h
   @brief Header for matcoal.c.
   @internal
 
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#ifndef MATCOAL_MATCOAL_H
#define MATCOAL_MATCOAL_H

#include "typedefs.h"
#include "pophist.h"

void     MatCoal_project(unsigned nSamples, double *x, double v, double errTol);
void     MatCoal_project_multi(unsigned nEpochs,unsigned nSamples,
							   double x[nEpochs][nSamples],
							   double tvec[nEpochs], PopHist *ph,
                               double errTol);
void MatCoal_integrate_epoch(unsigned nSamples,
                             double p0[nSamples],
                             double dt, double twoN, double errTol,
                             double p1[nSamples], double m[nSamples]);
void     MatCoal_integrate(unsigned nSamples, double m[nSamples], PopHist *ph,
                           double errTol);

#endif
