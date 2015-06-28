/**
 * @file model.c
 * @brief Functions for objects of type Model, a generic class whose
 * details are implemented elsewhere. The Model class generalizes the
 * task of calculating the expected value of an LD statistic from
 * parameter values involving recombination rate, c, mutation rate, u,
 * and the history of population size.
 *
 * Several of the models involve a vector of state variables, which is
 * iterated forward from generation to generation. This file
 * implements a function, odeEvolve, which approximates such systems as ODEs.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_odeiv2.h>
#include "misc.h"
#include "tokenizer.h"
#include "pophist.h"
#include "hill.h"
#include "strobeck.h"
#include "hayes.h"
#include "model.h"

pthread_mutex_t stderr_mutex = PTHREAD_MUTEX_INITIALIZER;

#if 0
pthread_mutex_t tmr_mutex = PTHREAD_MUTEX_INITIALIZER;
long        tmr_count;
clock_t     tmr_cputime;

#define TMR(E) tmr_begin=clock();\
    (E);\
    tmr_end=clock();\
    status=pthread_mutex_lock(&tmr_mutex);\
    if(status) ERR(status,"lock tmr_mutex");\
    ++tmr_count;\
    tmr_cputime += tmr_end - tmr_begin;\
    status=pthread_mutex_unlock(&tmr_mutex);\
    if(status) ERR(status,"unlock tmr_mutex");
#endif

struct ODE {
    const gsl_odeiv2_step_type *stepType;
    gsl_odeiv2_step *stepFun;
    gsl_odeiv2_control *ctrl;
    gsl_odeiv2_evolve *evolve;

    Model      *model;
    int         ydim;

    void       *data;
    struct dydt_params dydtPar;
};

static Model *Model_dup(const Model * model);

ODE        *ODE_new(const Model * model, double abstol, double reltol) {
    ODE        *ode = malloc(sizeof(ODE));

    checkmem(ode, __FILE__, __LINE__);

    ode->model = Model_dup(model);

    int         ydim = Model_stateDim(model);

    ode->ydim = ydim;

    ode->stepType = gsl_odeiv2_step_rk8pd;
    ode->stepFun = gsl_odeiv2_step_alloc(ode->stepType, ydim);
    ode->ctrl = gsl_odeiv2_control_y_new(abstol, reltol);
    ode->evolve = gsl_odeiv2_evolve_alloc(ydim);

    ode->data = Model_newState(model);
    memset(&ode->dydtPar, 0, sizeof(struct dydt_params));

    return ode;
}

void ODE_free(ODE * ode) {
    gsl_odeiv2_evolve_free(ode->evolve);
    gsl_odeiv2_control_free(ode->ctrl);
    gsl_odeiv2_step_free(ode->stepFun);

    Model_freeState(ode->model, ode->data);

    Model_free(ode->model);

    free(ode);
}

void       *ODE_state(ODE * ode) {
    return ode->data;
}

double ODE_stateVal(const ODE * ode, unsigned i) {
    return ode->model->stateVal(ode->data, i);
}

size_t ODE_stateDim(const ODE * ode) {
    return ode->ydim;
}

Model      *ODE_model(ODE * ode) {
    return ode->model;
}

struct dydt_params *ODE_dydtPar(ODE * ode) {
    return &ode->dydtPar;
}

double ODE_ld(ODE * ode, double c, double u, PopHist * ph) {
    assert(ode);
    return ode->model->ld(ode, c, u, ph, ode->model->twoNsmp);
}

/**
 * Allocate a new Model of type specified by the character string
 * "method".
 *
 * @param[in] method character string representing a method to be used
 * in predicting LD from population history. Legal values: "hill",
 * "strobeck", and "hayes".
 *
 * @param[in] twoNsmp haploid sample size, which is used to
 * incorporate sampling bias into predicted values. To turn this
 * feature off, use the value 0.
 *
 * @returns  a pointer to a newly-allocated Model, if successful; or
 * NULL, on failure.
 */
Model      *Model_alloc(const char *method, int twoNsmp) {
    Model      *model = NULL;
    char        buff[20];

    snprintf(buff, sizeof(buff), "%s", method);
    strlowercase(buff);
    if(strcmp(buff, "hill") == 0)
        model = Model_allocHill(twoNsmp);
    else if(strcmp(buff, "strobeck") == 0)
        model = Model_allocStrobeck(twoNsmp);
    else if(strcmp(buff, "hayes") == 0)
        model = Model_allocHayes(twoNsmp);
    else {
        free(model);
        return NULL;
    }

    return model;
}

/** Free an object of type Model */
void Model_free(Model * m) {
    free(m);
}

static Model *Model_dup(const Model * model) {
    Model      *m = malloc(sizeof(Model));

    memcpy(m, model, sizeof(Model));

    return m;
}

void       *Model_newState(const Model * m) {
    return m->newState();
}

void Model_freeState(Model * m, void *state) {
    m->freeState(state);
}

double Model_exactLD(const Model *model, ODE *ode, double c, double u,
                     PopHist *ph) {
    myassert(ode != NULL);
    return model->exactLD(c, u, ph, model->twoNsmp, ode->data);
}

/**
 * Calculate vector, ld, of values of sigdsq without using ODE
 * approximation.  
 *
 * @param[in] model An object of type Model, which determines
 *            which model to use in calculating LD.
 * @param[in,out] ode An object of type ODE, which is used only for
 *            its data field.
 * @param[in] nbins Size of arrays ld and c.
 * @param[out] ld is a vector of "nbins" doubles. On return, the i'th
 * entry will contain the value of sigma_d^2 implied by recombination
 * rate c[i], and by the population history in argument "ph".
 * @param[in] c Array of recombination rates.
 * @param[in] u Mutation rate.
 * @param[in] ph Population history.
 */
void Model_exactLDvec(const Model *model, ODE *ode, int nbins,
                        double ld[nbins], double c[nbins], double u,
                        PopHist *ph) {
    assert(model != NULL);
    assert(ode != NULL);
    int i;
    for(i=0; i < nbins; ++i)
        ld[i] = model->exactLD(c[i], u, ph, model->twoNsmp, ode->data);
}

/** Print state vector to file fp. */
void ODE_printState(const ODE * ode, FILE * fp) {
    unsigned    i, dim = ODE_stateDim(ode);

    fprintf(fp, "[");
    for(i = 0; i < dim; ++i) {
        fprintf(fp, "%lg", ODE_stateVal(ode, i));
        if(i + 1 < dim)
            fprintf(fp, ", ");
    }
    fprintf(fp, "]");
}

/**
 * Evolve through all Epochs of population history, beginning with
 * initial state given in y. Calculation uses ODE approximation. The
 * initial value of y should be set before calling ODE_evolve.
 *
 * @param[in,out] ode an object of type ODE.
 * @param[in,out] y Vector of state variables.
 * @param[in] odeStepSize Controls the size of step taken by the minimizer.
 * @param[in] ph Describes the population's history. If ph contains
 *            just one Epoch, y gets its equilibrium value.
 * @param[in] verbose Verbosity.
 * @param[in] c Recombination rate; overrides values in ph.
 * @param[in] u Mutation rate.
 * @returns Returns 0 on success, 1 if function didn't run because
 *            dydt was not provided. 
 */
int ODE_evolve(ODE * ode,
               double *y,
               unsigned ydim,
               double c,
               double u,
               PopHist * ph,
               int dydt(double t_notused, const double y1[], double f[],
                        void *params), int verbose) {

    int         i;
    double      t = 0.0;
    double      stepSize;
    gsl_odeiv2_system sys = { dydt, NULL, ydim, &ode->dydtPar };

    ode->dydtPar.c = c;
    ode->dydtPar.u = u;

#ifdef DEBUG
    assertFiniteArray(y, ydim, __FILE__, __LINE__);
#endif

    gsl_odeiv2_step_reset(ode->stepFun);
    gsl_odeiv2_evolve_reset(ode->evolve);

    for(i = PopHist_nepoch(ph) - 2; i >= 0; --i) {
        long unsigned nsteps = 0;
        double      t1 = PopHist_duration(ph, i);

        assert(isfinite(t1));
        ode->dydtPar.duration = t1;
        ode->dydtPar.twoN = PopHist_twoN(ph, i);
        stepSize = 0.05 * t1;

        while(t < t1) {
#ifdef DEBUG
            assertFiniteArray(y, ydim, __FILE__, __LINE__);
#endif
            int         status;

#if 0
            status = pthread_mutex_lock(&stderr_mutex);
            if(status)
                ERR(status, "lock stderr_mutex");
            fprintf(stderr, "t=%lg t1=%lg stepSize=%lg\n", t, t1, stepSize);
            status = pthread_mutex_unlock(&stderr_mutex);
            if(status)
                ERR(status, "lock stderr_mutex");
#endif

            /*
             * Following line dumps core because, at some point, t >
             * t1, but ode->stepsize > 0. This shouldn't happen
             * because of the condition in the while loop above. Some
             * other thread must be writing to t or to t1.
             */
            status = gsl_odeiv2_evolve_apply(ode->evolve,
                                             ode->ctrl,
                                             ode->stepFun,
                                             &sys, &t, t1, &stepSize, y);

#ifdef DEBUG
            assertFiniteArray(y, ydim, __FILE__, __LINE__);
#endif
            ++nsteps;           /* debug */

            if(verbose) {
                int         j;

                printf("ODE_evolve step %lu: curr state =", nsteps);
                for(j = 0; j < ydim; ++j)
                    printf(" %lg", y[j]);
                printf(" twoN=%lg t=%lg", ode->dydtPar.twoN,
                       ode->dydtPar.duration);
                putchar('\n');
            }

            if(status != GSL_SUCCESS)
                break;
        }

#ifdef DEBUG
        assertFiniteArray(y, ydim, __FILE__, __LINE__);
#endif
    }

    return 0;
}

/**
 * On entry, c is a list of recombination rates, ph the population
 * history, and odeStepSize is the initial step size.
 *
 * @param[in,out] ode An object of type ODE.
 * @param[in] c,u Rates of recombination and mutation.
 * @param[in] nbins Size of arrays sigdsq and c.
 * @param[in] ph Population history.
 * @param[out] ld is a vector of "nbins" doubles. On return, the i'th
 * entry will contain the value of sigma_d^2 implied by recombination
 * rate c[i], and by the population history in argument "ph".
 */
void ODE_ldVec(ODE * ode, double *ld, int nbins, const double *c, double u,
               PopHist * ph) {
    int         i;

    for(i = 0; i < nbins; ++i) {
        ld[i] = ode->model->ld(ode, c[i], u, ph, ode->model->twoNsmp);
    }

    return;
}

double ODE_ldEq(ODE * ode, double c, double u, PopHist * ph,
                unsigned whichEpoch) {
    assert(ode);
    return ode->model->ldEq(c, u, ph, whichEpoch, ode->model->twoNsmp,
                            ode->data);
}

/**
 * This function returns 1 vector:
 * eq[i] is sigma_d^2 for a population at equilibrium with parameters
 *       as given in the specified PopHist Epoch, and recombination rate as
 *       specified by c[i].
 *
 * On entry, c is a list of recombination rates and ep the population
 * history Epoch.
 *
 * @param[in,out] model An object of type Model, which specifies
 * the method to be used in calculating expected LD.
 */
void ODE_ldVecEq(ODE * ode, double *eq, int nbins, const double *c,
                 double u, PopHist * ph, int whichEpoch) {
    int         i;

    for(i = 0; i < nbins; ++i)
        eq[i] = ode->model->ldEq(c[i], u, ph, whichEpoch,
                                 ode->model->twoNsmp, ode->data);

    return;
}

/**
 * Allocate a new ModelList.
 *
 * @param[in] s0 A comma-separated list of Model types, such as
 * "Hill,Strobeck".
 */
ModelList  *ModelList_alloc(const char *s0, int twoNsmp) {
    char        buff[100];
    Tokenizer  *tkz = Tokenizer_new(MAXMODELS);
    unsigned    ntokens, token;
    ModelList  *ml = malloc(sizeof(ModelList));

    checkmem(ml, __FILE__, __LINE__);
    ml->nModels = 0;

    snprintf(buff, sizeof(buff), "%s", s0);

    Tokenizer_split(tkz, buff, ",");
    ntokens = Tokenizer_strip(tkz, " \t\r\"");

    for(token = 0; token < ntokens; ++token)
        ModelList_addModel(ml, Tokenizer_token(tkz, token), twoNsmp);

    Tokenizer_free(tkz);
    return ml;
}

/**
 * Add a method to an ModelList. Return 1 on success, 0 on failure.
 *
 * @param[out] ml The ModelList.
 * @param[in] modelName A character string such as "hill", or
 * "strobeck".
 * @param[in] twoNsmp The number of individuals sampled.
 */
unsigned ModelList_addModel(ModelList * ml, const char *modelName,
                            int twoNsmp) {

    Model      *model = Model_alloc(modelName, twoNsmp);

    if(model == NULL)
        return 0;
    ml->m[ml->nModels++] = model;
    return 1;
}

void ModelList_free(ModelList * el) {
    unsigned    i;

    for(i = 0; i < el->nModels; ++i)
        Model_free(el->m[i]);
    free(el);
    return;
}

