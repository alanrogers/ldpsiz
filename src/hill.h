/**
 * @file hill.h
 * @author Alan R. Rogers
 * @brief Header for hill.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_HILL
#  define LDPSIZ_HILL

#  include "typedefs.h"
#  include "pophist.h"

void       *HillData_new(void);
double      HillData_stateVal(void *vdata, unsigned i);
void        HillData_free(void *p);

int         Hill_geteq(double x[], double twoN, double c, double u);
int         Hill_dydt(double t, const double y[], double f[], void *params);
void        Hill_approxeq(double x[], double twoN, double c, double u);
double      Hill_sigdsq(ODE * ode, double c, double u, PopHist * ph,
                        int twoNsmp);
double      Hill_sigdsqEq(double c, double u, PopHist * ph,
                          unsigned whichEpoch, int twoNsmp, void *data);
size_t      Hill_stateDim(void);
const char *Hill_stateLbl(unsigned i);

Model      *Model_allocHill(int twoNsmp);
int         Hill_evolveDiscrete(double y[], PopHist * ph, double c, double u);
double      Hill_get_sigdsq(double y[], int twoNsmp);

#  ifndef NDEBUG
int         Hill_test_getDRM(double twoN, double c, double u);
int         Hill_test_onestep(double *x1, int dim, double twoN, double c, double u);
int         onestep_slow(const double x[], double y[], double twoN, double c,
						 double u);
#  endif

#endif
