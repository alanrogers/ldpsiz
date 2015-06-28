/**
 * @file model.h
 * @author Alan R. Rogers
 * @brief Header for model.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_MODEL
#define LDPSIZ_MODEL

#include <assert.h>
#include "typedefs.h"
#include "pophist.h"
#include "misc.h"

#define MODEL_LBL_SIZE 10
#define MAXMODELS 3

#if 1
   /* Hayes method assumes E[rsq] = 1/(4*N*c + 2) */
#  define HAYES_MUTATION_ADJUSTMENT
#else
   /* Hayes method assumes E[rsq] = 1/(4*N*c + 1) */
#  undef HAYES_MUTATION_ADJUSTMENT
#endif

/**
 * An object of type Model is shared among threads and among tasks. It
 * therefore cannot have any per-thread or per-task data.
 */
struct Model {
    char        lbl[MODEL_LBL_SIZE];
    int         twoNsmp;        /* number of gene copies in sample */

    /* calc expected LD */
    double      (*ld) (ODE * ode, double c, double u, PopHist * ph,
                       int twoNsmp);

    /* calc expected equilibrium LD */
    double      (*ldEq) (double c, double u, PopHist * ph,
                         unsigned whichEpoch, int twoNsmp, void *data);

    /* calc LD without using ODE */
    double      (*exactLD)(double c, double u, PopHist *ph, int twoNsmp,
                           void *state);

    /* return dimension of state vector */
                size_t(*stateDim) (void);

    /* return i'th state label */
    const char *(*stateLbl) (unsigned i);

    /* return i'th state value */
    double      (*stateVal) (void *vdata, unsigned i);

    /* allocate state variables */
    void       *(*newState) (void);

    /* free state variables */
    void        (*freeState) (void *vdata);
};

/** A list of models */
struct ModelList {
    unsigned    nModels;
    Model      *m[MAXMODELS];
};

double      Model_exactLD(const Model *model, ODE *ode, double c, double u,
                          PopHist *ph);
void Model_exactLDvec(const Model *model, ODE *ode, int nbins,
                      double ld[nbins], double c[nbins], double u,
                      PopHist *ph);
Model      *Model_alloc(const char *method, int twoNsmp);
void        Model_free(Model * m);
static inline const char *Model_lbl(const Model * model);
void       *Model_newState(const Model * m);
void        Model_freeState(Model * m, void *state);
static inline size_t Model_stateDim(const Model * model);
static inline const char *Model_stateLbl(const Model * model, unsigned i);
unsigned    ModelList_addModel(ModelList * el, const char *modelName,
                               int twoNsmp);
ModelList  *ModelList_alloc(const char *s0, int twoNsmp);
void        ModelList_free(ModelList * el);
static inline Model *ModelList_model(ModelList * el, unsigned ndx);
static inline unsigned ModelList_size(const ModelList * el);

/* inline functions */
static inline const char *Model_lbl(const Model * model) {
    assert(model);
    return model->lbl;
}

static inline const char *Model_stateLbl(const Model * model, unsigned i) {
    if(model->stateLbl)
        return model->stateLbl(i);
    else
        eprintf("ERR@:%s:%d: There are no state labels for model \"%s\"",
                __FILE__, __LINE__, model->lbl);
    return NULL;
}

static inline size_t Model_stateDim(const Model * model) {
    assert(model);
    return model->stateDim();
}

static inline unsigned ModelList_size(const ModelList * el) {
    assert(el);
    return el->nModels;
}

static inline Model *ModelList_model(ModelList * ml, unsigned ndx) {
    assert(ml);
    assert(ndx < ml->nModels);
    return ml->m[ndx];
}

ODE        *ODE_new(const Model * model, double abstol, double reltol);
void        ODE_free(ODE * ode);
void       *ODE_state(ODE * ode);
double      ODE_stateVal(const ODE * ode, unsigned i);
size_t      ODE_stateDim(const ODE * ode);
Model      *ODE_model(ODE * ode);
struct dydt_params *ODE_dydtPar(ODE * ode);
void        ODE_printState(const ODE * ode, FILE * fp);
void        ODE_ldVec(ODE * ode, double *ld, int nbins, const double *c,
                      double u, PopHist * ph);
double      ODE_ldEq(ODE * ode, double c, double u, PopHist * ph,
                     unsigned whichEpoch);
void        ODE_ldVecEq(ODE * ode, double *eq, int nbins, const double *c,
                        double u, PopHist * ph, int whichEpoch);
double      ODE_ld(ODE * ode, double c, double u, PopHist * ph);
int         ODE_evolve(ODE * ode,
                       double *y,
                       unsigned ydim,
                       double c,
                       double u,
                       PopHist * ph,
                       int dydt(double t_notused, const double y1[],
                                double f[], void *params), int verbose);

#endif
