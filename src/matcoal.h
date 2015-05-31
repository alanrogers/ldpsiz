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

void MatCoal_project(unsigned nSamples, long double *x,
					 long double v, long double betavec[nSamples],
					 double errTol);
void     MatCoal_project_multi(unsigned nEpochs,unsigned nSamples,
							   double x[nEpochs][nSamples],
							   double tvec[nEpochs], PopHist *ph,
                               double errTol);
void MatCoal_integrate_epoch(unsigned nSamples,
                             long double p0[nSamples],
                             double dt, double twoN, double errTol,
                             long double p1[nSamples],
							 double m[nSamples],
                             long double betavec[nSamples]);
void     MatCoal_integrate(unsigned nSamples, double m[nSamples], PopHist *ph,
                           double errTol);
static inline double MatCoal_beta(unsigned i);

/**
 * Used to construct the entries of the transition rate matrix, B.
 * B[i][i] = -beta(i) and B[i][i+1] = beta(i+1).  (All the other
 * entries of B are zero.)  In this code, index i corresponds to the
 * coalescent interval containing i+1 distinct lineages.  Thus, the
 * code below corresponds to j*(j-1)/2, where j is the number of
 * lineages in the coalescent interval.
 */
static inline double MatCoal_beta(unsigned i) { return (i*(i+1))/2; }

#endif
