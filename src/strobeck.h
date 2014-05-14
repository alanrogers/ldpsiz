/**
 * @file strobeck.h
 * @author Alan R. Rogers
 * @brief Header for strobeck.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_STROBECK
#define LDPSIZ_STROBECK

#include "typedefs.h"
#include "pophist.h"

void       *StrobeckData_new(void);
void        StrobeckData_free(void *p);
double      StrobeckData_stateVal(void *vdata, unsigned i);

int         Strobeck_geteq(double x[], double twoN, double c, double u);
int         Strobeck_dydt(double t, const double y[], double f[],
                          void *params);
double      Strobeck_get_sigdsq(double y[], int n);
int         Strobeck_evolveDiscrete(double y[], PopHist * ph, double c,
                                    double u);
double      Strobeck_sigdsq(ODE * ode, double c, double u, PopHist * ph,
                            int twoNsmp);
double      Strobeck_sigdsqEq(double c, double u, PopHist * ph,
                              unsigned whichEpoch, int twoNsmp, void *data);
size_t      Strobeck_stateDim(void);
const char *Strobeck_stateLbl(unsigned i);
Model      *Model_allocStrobeck(int twoNsmp);

#endif
