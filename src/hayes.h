/**
 * @file fithayes.h
 * @author Alan R. Rogers
 * @brief Header for hayes.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_HAYES_H
#define LDPSIZ_HAYES_H

#include "pophist.h"

double      Hayes_rsqNoODE(double c, double u, PopHist *ph, int twoNsmp,
                           void *notused);
double      Hayes_rsq(ODE * unused0, double c, double u, PopHist * ph,
                      int twoNsmp);
double      Hayes_rsqEq(double c, double u, PopHist * ph, unsigned whichEpoch,
                        int twoNsmp, void *notused);
size_t      Hayes_stateDim(void);
const char *Hayes_stateLbl(unsigned notused);
void       *HayesData_new(void);
double      HayesData_stateVal(void *notused, unsigned i);
void        HayesData_free(void *notused);

Model      *Model_allocHayes(int twoNsmp);
#endif
