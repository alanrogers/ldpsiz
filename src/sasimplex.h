/* 
 * Copyright (C) 2014 Alan Rogers
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
#ifndef __GSL_SASIMPLEX_H__
#define __GSL_SASIMPLEX_H__

#include <gsl/gsl_types.h>
#include <gsl/gsl_multimin.h>

GSL_VAR const gsl_multimin_fminimizer_type *gsl_multimin_fminimizer_sasimplex;

void        sasimplex_random_seed(gsl_multimin_fminimizer * minimizer,
                                  unsigned seed);
void        sasimplex_set_temp(gsl_multimin_fminimizer * minimizer,
                               double temperature);
int         sasimplex_set_bounds(gsl_multimin_fminimizer * minimizer,
                                 const gsl_vector *lbound,
                                 const gsl_vector *ubound,
                                 const gsl_vector *step_size);
int         sasimplex_randomize_state(gsl_multimin_fminimizer * minimizer,
                                      int rotate, gsl_vector * lo,
                                      gsl_vector * hi,
                                      const gsl_vector * step_size);
double      sasimplex_vertical_scale(gsl_multimin_fminimizer * minimizer);
int         sasimplex_n_iterations(gsl_multimin_fminimizer * minimizer,
                                   double *size,
                                   double tol_fval,
                                   double tol_size,
                                   int nItr, double temperature, int verbose);
int         sasimplex_converged(gsl_multimin_fminimizer * minimizer,
                                double *size, double tol_fval,
                                double tol_size);
void        sasimplex_print(gsl_multimin_fminimizer * minimizer);

#endif                       /* __GSL_SASIMPLEX_H__ */
