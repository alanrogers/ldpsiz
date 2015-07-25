/**
 * @file pophist.h
 * @author Alan R. Rogers
 * @brief Header for pophist.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_POPHIST
#define LDPSIZ_POPHIST

#define PERTURB_GAUSSIAN
#define PAR_PER_EPOCH 2 ///< The number of parameters in each Epoch

/* #define PERTURB_TDIST */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stddef.h>
#include "typedefs.h"

/**
 * PopHist represents the history of a single population.
 * It is a structure allocated on a single block of memory, using the
 * "struct hack" of C programming. This makes it possible to avoid
 * pointers within PopHist, so that these objects can be copied using
 * memcpy. (Otherwise, you end up with pointers in one object pointing
 * to memory in another.)
 */
struct PopHist {
    unsigned    nepoch;    /**< number of epochs */
    size_t      size;      /**< size of memory block allocated */
    double      p[PAR_PER_EPOCH];   /**< array of parameter values */
};

EpochLink  *EpochLink_new(EpochLink * head, double t, double n);
EpochLink  *EpochLink_dup(EpochLink * head);
void        EpochLink_free(EpochLink * link);
void        EpochLink_print(int ndx, const EpochLink * link, FILE * ofp);
int         EpochLink_nlinks(const EpochLink * head);
const EpochLink *EpochLink_next(const EpochLink * link);
double      EpochLink_duration(const EpochLink * link);
double      EpochLink_twoN(const EpochLink * link);
#ifndef NDEBUG
void        EpochLink_test(void);
#endif

void        PopHist_init(PopHist * ph, unsigned nepoch, size_t size);
PopHist    *PopHist_newEmpty(unsigned nepoch);
size_t      PopHist_calcSize(unsigned nepoch);
double      PopHist_age(const PopHist *ph, int i);
int         PopHist_findEpoch(const PopHist *ph, double age);
double      PopHist_twoNinv(const PopHist * ph, int i);
void        PopHist_setDuration(PopHist * ph, int i, double duration);
void        PopHist_setTwoN(PopHist * ph, int i, double twoN);
void        PopHist_setTwoNinv(PopHist * ph, int i, double twoNinv);
PopHist    *PopHist_fromEpochLink(const EpochLink * head);
unsigned    PopHist_nepoch(const PopHist * ph);
void        PopHist_print_comment(const PopHist * ph, const char *comstr,
                                  FILE * outfile);
void        PopHist_print(const PopHist * ph, FILE * outfile);
void        PopHist_free(PopHist * ph);
unsigned    PopHist_nParams(const PopHist * ph);
int         PopHist_paramName(const PopHist * ph, char *buff, int bufflen,
                              int ndx);
double      PopHist_paramValue(const PopHist * ph, int ndx);
void        PopHist_to_vector(gsl_vector * v, const PopHist * ph);
void        PopHist_to_C_array(int dim, double v[dim], const PopHist * ph);
void        PopHist_to_C_array_invert_N(double *v, int dim,
                                        const PopHist * ph);
void        C_array_to_PopHist(PopHist * ph, int dim, const double v[dim]);
void        PopHist_setSimplexScale(gsl_vector * scale, const PopHist * ph);
void        PopHist_setAllTwoNinv(double *x, int dim, double value);
void        PopHist_setAllDuration(double *x, int dim, double value);
void        vector_to_PopHist(PopHist * ph, const gsl_vector * v);
PopHist    *PopHist_dup(const PopHist * ph);
void        PopHist_copy(PopHist * dest, const PopHist * src);
PopHist    *PopHist_perturb(PopHist * ph, const PopHist * model, double dt,
                            double dN, gsl_rng * rng);
PopHist    *PopHist_perturbInPlace(PopHist * ph, double dt, double dN,
                                   gsl_rng * rng);
void        PopHist_sanityCheck(PopHist * ph, const char *file, int line);
double      PopHist_distance(PopHist * ph1, PopHist * ph2);
size_t      PopHist_size(const PopHist * ph);
int         PopHist_calc_nParams(int nepoch);

static inline double PopHist_duration(const PopHist * ph, int i);
static inline double PopHist_twoN(const PopHist * ph, int i);

/// Return duration of Epoch i.
static inline double PopHist_duration(const PopHist * ph, int i) {
    return ph->p[1 + i * PAR_PER_EPOCH];
}

/// Return the haploid population size, 2N, during epoch i.
static inline double PopHist_twoN(const PopHist * ph, int i) {
    return 1.0 / ph->p[i * PAR_PER_EPOCH];
}
#endif
