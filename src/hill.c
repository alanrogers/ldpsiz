/**
 * @file hill.c
 * @author Alan R. Rogers
 * @brief Method of Hill, William G. 1975. TPB 8:117-126.
 *
 * The state variable is as defined on p 120 of Hill (1975). Hill's x
 * corresponds to entries 0, 1, and 2 of my vector below. Hill also
 * needs the expected heterozygosity (see his eqn 7), so I have added
 * that on as element 3 of the state vector.
 *
 * y[0] = E[H_A*H_B]
 * y[1] = 4 E[sum_{ij} a_i b_j D_{ij}]
 * y[2] = 2 E[sum_{ij} D_{ij}^2
 * y[3] = E[H_A] = E[H_B]
 *
 * Here, a_i and b_j are the frequencies of the i'th allele at locus A
 * and the j'th allele at locus B. D_{ij} is the coefficient of
 * linkage disequilibrium between these alleles, H_A = 1 - sum a_i^2
 * is the heterozygosity at locus A and H_B is the corresponding value
 * for locus B.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "misc.h"
#include "pophist.h"
#include "model.h"
#include "hill.h"

#define DIM 4

#undef DEBUG

#undef DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

static const char *stateLbl[DIM] = { "Ha*Hb", "4S[abD]", "2S[D^2]", "H" };

struct HillData {
    unsigned    ydim;
    double      y[DIM];
};

static int  onestep(const double x[], double y[], double twoN, double c,
                    double u);
static int  getDRM(double DRM[][DIM], double twoN, const double c, double u);

#ifndef NDEBUG
static int  getDRM_slow(double DRM[][DIM], double twoN, const double c,
						double u);
#endif

int         iterate(double y[], int t, double twoN, double c, double u);
double      Hill_exactLD(double c, double u, PopHist *ph, int twoNsmp,
                         void *data);

void       *HillData_new(void) {
    HillData   *hd = malloc(sizeof(HillData));

    checkmem(hd, __FILE__, __LINE__);
    memset(hd, 0, sizeof(HillData));
    hd->y[0] = 1.0;
    hd->y[1] = 2.0;
    hd->y[2] = 3.0;
    hd->y[3] = 4.0;
    hd->ydim = DIM;
    return hd;
}

/** i'th state value */
double HillData_stateVal(void *vdata, unsigned i) {
    HillData   *data = vdata;

    return data->y[i];
}

/** free state variables */
void HillData_free(void *vdata) {
    HillData   *data = vdata;

    free(data);
}

Model      *Model_allocHill(int twoNsmp) {
    Model      *model = malloc(sizeof(Model));

    checkmem(model, __FILE__, __LINE__);

    snprintf(model->lbl, sizeof(model->lbl), "Hill");
    model->twoNsmp = twoNsmp;

    model->ld = Hill_sigdsq;
    model->exactLD = Hill_exactLD;
    model->ldEq = Hill_sigdsqEq;
    model->stateDim = Hill_stateDim;
    model->stateLbl = Hill_stateLbl;
    model->stateVal = HillData_stateVal;
    model->newState = HillData_new;
    model->freeState = HillData_free;

    return model;
};

/* dimension of state vector */
size_t Hill_stateDim(void) {
    return DIM;
}

const char *Hill_stateLbl(unsigned i) {
    return stateLbl[i];
}

/*
 * Uses definitions on Hill's p 119.
 */
static int getDRM(double DRM[][DIM], double twoN, const double c, double u) {
    double      N = 0.5 * twoN; /* number of diploids individuals */
    double      m, m2, m4, one_minus_c;

    /* Construct RM vector (recombination and mutation) */
    m = (1.0 - u);
    assert(m > 0.0);
    m2 = m * m;
    m4 = m2 * m2;
    one_minus_c = 1.0 - c;
    assert(one_minus_c > 0.0);
    double      RM[DIM] = { m4,
        one_minus_c * m4,
        one_minus_c * one_minus_c * m4,
        m2
    };

    /*
     * Construct drift matrix, DRM = (I - P/N)*Diag(RM)
     */
    DRM[0][0] = RM[0] * (1.0 - 2.0 / twoN);
    DRM[0][1] = RM[1] / twoN;
    DRM[0][2] = 0.0;
    DRM[0][3] = 0.0;

    DRM[1][0] = 0.0;
    DRM[1][1] = RM[1] * (1.0 - 5.0 / twoN);
    DRM[1][2] = RM[2] / N;
    DRM[1][3] = 0.0;

    DRM[2][0] = RM[0] / N;
    DRM[2][1] = RM[1] / N;
    DRM[2][2] = RM[2] * (1.0 - 3.0 / twoN);
    DRM[2][3] = 0.0;

    DRM[3][0] = 0.0;
    DRM[3][1] = 0.0;
    DRM[3][2] = 0.0;
    DRM[3][3] = RM[3] * (1.0 - 1.0 / twoN);

    return 0;
}

/*
 * Find equilibrium vector, x.  Hill's Eqn 11 implies that x is the
 * solution of
 *
 *   (I - DRM)*x = w
 *
 * for known DRM and w.  This linear system is solved here using Gaussian
 * elimination.
 */
int Hill_geteq(double x[], double twoN, double c, double u) {
    unsigned    i;
    int         s;
    double      DRM[DIM][DIM];
    double      U = 0.5 * twoN * u;
    double      z;

    /* Equilibrium w is (Hill's Eqn 8). */
    double      w_data[DIM] = { 16.0 * U * u / (4.0 * U + 1.0),
        0.0, 0.0, 2.0 * u
    };
    gsl_vector_view w = gsl_vector_view_array(w_data, DIM);

    gsl_vector_view xvec = gsl_vector_view_array(x, DIM);

    gsl_matrix_view A = gsl_matrix_view_array(&DRM[0][0], DIM, DIM);
    gsl_permutation *p = gsl_permutation_alloc(DIM);

    getDRM(DRM, twoN, c, u);

    /* Form A = I-DRM */
    gsl_matrix_scale(&A.matrix, -1.0);
    for(i = 0; i < DIM; ++i) {
        z = gsl_matrix_get(&A.matrix, i, i);
        gsl_matrix_set(&A.matrix, i, i, 1.0 + z);
    }

#if 0
    int         status;

    status = pthread_mutex_lock(&outputLock);
    if(status)
        ERR(status, "lock output");

    fprintf(stderr, "%s:%s:%d: %s", __FILE__, __func__, __LINE__,
            " GSL Linear system (I - DRM) x = w:\n");
    for(i = 0; i < DIM; ++i) {
        putc('[', stderr);
        for(unsigned j = 0; j < DIM; ++j)
            fprintf(stderr, " %8.5f", gsl_matrix_get(&A.matrix, i, j));
        fprintf(stderr, "] [x_%d]  [%8.5g]\n", i,
                gsl_vector_get(&w.vector, i));
    }

    status = pthread_mutex_unlock(&outputLock);
    if(status)
        ERR(status, "unlock output");
#endif

    /* Solve linear system using LU decomposition */
    gsl_permutation_init(p);
    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_solve(&A.matrix, p, &w.vector, &xvec.vector);

#if 0
    status = pthread_mutex_lock(&outputLock);
    if(status)
        ERR(status, "lock output");

    fprintf(stderr, "%s:%s:%d: Solution vector:",
            __FILE__, __func__, __LINE__);
    for(i = 0; i < DIM; ++i)
        fprintf(stderr, " %lg", x[i]);
    putc('\n', stderr);

    status = pthread_mutex_unlock(&outputLock);
    if(status)
        ERR(status, "unlock output");

#endif

    gsl_permutation_free(p);

    return 0;
}

/* Approximate equilibrium.  Hill's equations 7 and 10 */
void Hill_approxeq(double x[], double twoN, double c, double u) {
    double      U = 0.5 * twoN * u;
    double      C = 0.5 * twoN * c;
    double      V = U;
    double      UpV = U + V;

    double      K = 8 * U * V * (1 / (1 + 4 * U) + 1 / (1 + 4 * V)) /
        (9 + 26 * C + 54 * UpV + 8 * C * C + 76 * C * UpV
         + 80 * UpV * UpV + 16 * C * C * UpV + 48 * C * UpV * UpV
         + 32 * UpV * UpV * UpV);

    x[0] = 11 + 26 * C + 32 * UpV + 8 * C * C + 24 * C * UpV + 16 * UpV * UpV;
    x[1] = 4;
    x[2] = 10 + 4 * C + 8 * UpV;

    x[0] *= K;
    x[1] *= K;
    x[2] *= K;

    x[3] = 4 * U / (4 * U + 1); /* Hill's eqn 7 */

}

/** Iterate the difference equation across t generations. */
int iterate(double y[], int t, double twoN, double c, double u) {
    double      y2[DIM];

    while(t-- > 0) {
        onestep(y, y2, twoN, c, u);
        memcpy(y, y2, sizeof(y2));
    }
    return 0;
}

/**
 * Move the state vector forward one generation, using Hill's
 * recurrence (Eqns 5 and 8).
 *
 * @param[input] x is the initial state vector.
 * @param[output] y is the new state vector.
 */
static int onestep(const double x[], double y[], double twoN, double c,
                   double u) {
    /*double U = twoN*u; */
    double      DRM[DIM][DIM];

#ifdef DEBUG
    assertFiniteArray(x, DIM, __FILE__, __LINE__);
#endif

    /* DRM is matrix of effects of drift, recombination and mutation */
    getDRM(DRM, twoN, c, u);

#ifdef DEBUG
    if(!matIsFinite(DIM, DRM)) {
        printf("onestep got nonfinite DRM matrix.\n");
        printf("  on entry: [%lg, %lg, %lg, %lg]\n", x[0], x[1], x[2], x[3]);
        printsqrmat("nonfinite DRM", DIM, DRM);
        exit(1);
    }
#endif

    /* begin with increment due to mutation (Hill's Eqn 8). */
    y[0] = 4.0 * u * x[3];
    y[1] = y[2] = 0.0;
    y[3] = 2.0 * u;

    /*
     * Add the matrix-vector product, DRM*y (Hill's eqn 5).
     * Unrolled for speed. This had a huge effect.
     */
    y[0] +=
        DRM[0][0] * x[0] + DRM[0][1] * x[1] + DRM[0][2] * x[2] +
        DRM[0][3] * x[3];
    y[1] +=
        DRM[1][0] * x[0] + DRM[1][1] * x[1] + DRM[1][2] * x[2] +
        DRM[1][3] * x[3];
    y[2] +=
        DRM[2][0] * x[0] + DRM[2][1] * x[1] + DRM[2][2] * x[2] +
        DRM[2][3] * x[3];
    y[3] +=
        DRM[3][0] * x[0] + DRM[3][1] * x[1] + DRM[3][2] * x[2] +
        DRM[3][3] * x[3];

#ifdef DEBUG
    assertFiniteArray(y, DIM, __FILE__, __LINE__);
#endif

    return 0;
}

/**
 * params points to a structure of type dydt_params.
 *
 * On return, f[] contains the derivatives of y1[], as approximated by
 * difference equations.
 */
int Hill_dydt(double t_notused, const double x[], double f[], void *params) {
    double      y[DIM];
    struct dydt_params *par = (struct dydt_params *) params;

#ifdef DEBUG
    assertFiniteArray(x, DIM, __FILE__, __LINE__);
#endif

    onestep(x, y, par->twoN, par->c, par->u);

#ifdef DEBUG
    assertFiniteArray(y, DIM, __FILE__, __LINE__);
#endif

    /*
     * f is the difference between new and old state vectors.
     * Unrolled for speed.
     */
    f[0] = y[0] - x[0];
    f[1] = y[1] - x[1];
    f[2] = y[2] - x[2];
    f[3] = y[3] - x[3];

#ifdef DEBUG
    assertFiniteArray(f, DIM, __FILE__, __LINE__);
#endif

    return GSL_SUCCESS;
}

#ifndef NDEBUG
int Hill_test_getDRM(double twoN, double c, double u) {
    double      DRM1[DIM][DIM], DRM2[DIM][DIM];
    int i, j;

    getDRM(DRM1, twoN, c, u);
    getDRM_slow(DRM2, twoN, c, u);

    for(i = 0; i < DIM; ++i) {
        for(j = 0; j < DIM; ++j) {
            if(fabs(DRM1[i][j] - DRM2[i][j]) > DBL_EPSILON)
                return 1;
        }
    }
    return 0;
}

int Hill_test_onestep(double *x1, int dim, double twoN, double c, double u) {
    double      x2[DIM], yy1[DIM], y2[DIM];
    int i;
               
    assert(dim == DIM);

    memcpy(yy1, x1, DIM * sizeof(x1[0]));
    onestep_slow(x1, x2, twoN, c, u);
    onestep(yy1, y2, twoN, c, u);
    for(i = 0; i < DIM; ++i) {
        if(fabs(x2[i] - y2[i]) > DBL_EPSILON) 
            return 1;
    }
    return 0;
}

#endif


#if 0
static int  check_equilibrium(double y[], double tol, double c, double u,
                              double twoN, int verbose);

/* Compare equilibrium y with approximate formula */
static int check_equilibrium(double y[], double tol, double c, double u,
                             double twoN, int verbose) {
#if 0
    const char *lbl[] = { "H_A*H_B", "4SumpqD", "2SumDsq",
        "H"
    };
#endif
    int         rval = 0, nsteps = 1000;
    double      y2[DIM];
    double      relerr_it, relerr_approx;

    /*
     * Test by iterating difference equation.
     */
    memcpy(y2, y, sizeof(y2));
    iterate(y2, nsteps, twoN, c, u);
    relerr_it = getreldiff(DIM, y, y2, verbose);
    if(relerr_it > tol) {
        printf("check_equilibrium: %d-step test FAILED. Relerr=%lg\n",
               nsteps, relerr_it);
        rval = 1;
    }

    /*
     * Test against Hill's approximate formula
     */
    Hill_approxeq(y2, twoN, c, u);
    relerr_approx = getreldiff(DIM, y, y2, verbose);
    if(relerr_approx > tol) {
        printf("check_equilibrium: approx test FAILED. Relerr=%lg\n",
               relerr_approx);
        rval = 1;
    }
    return rval;
}
#endif

/**
 * Set state vector to equilibrium for Epoch "whichEpoch" and then use it
 * to calculate sigmdsq.
 */
double Hill_sigdsqEq(double c, double u, PopHist * ph, unsigned whichEpoch,
                     int twoNsmp, void *vdata) {

    HillData   *data = (HillData *) vdata;

    /* y is equilibrium for appropriate epoch */
    Hill_geteq(data->y, PopHist_twoN(ph, whichEpoch), c, u);

    return Hill_get_sigdsq(data->y, twoNsmp);
}

/**
 * Set state vector to initial equilibrium, then evolve it through the
 * population history, and then use this state vector to calculate sigmdsq.
 */
double Hill_sigdsq(ODE * ode, double c, double u, PopHist * ph, int twoNsmp) {
    HillData   *data = (HillData *) ODE_state(ode);

#ifdef DEBUG
    int         status, i;
    char        buff[20];

    for(i = 0; i < PopHist_nParams(ph); ++i) {
        PopHist_paramName(ph, buff, sizeof(buff), i);
    }

#endif

    /* Initial y is equilibrium for earliest epoch */
    unsigned    nepoch = PopHist_nepoch(ph);

    Hill_geteq(data->y, PopHist_twoN(ph, nepoch - 1), c, u);

    DPRINTF(("%s:%s:%d: y=[%lg, %lg, %lg, %lg]\n",
             __FILE__, __func__, __LINE__,
             data->y[0], data->y[1], data->y[2], data->y[3]));

    /* evolve through population history */
    ODE_evolve(ode, data->y, data->ydim, c, u, ph, Hill_dydt, !VERBOSE);

    DPRINTF(("%s:%s:%d: y=[%lg, %lg, %lg, %lg]\n",
             __FILE__, __func__, __LINE__,
             data->y[0], data->y[1], data->y[2], data->y[3]));

    double sigdsq = Hill_get_sigdsq(data->y, twoNsmp);
    return sigdsq;
}

/**
 * Calculate sigma_d^2 from state vector. On page 124, just after eqn 14,
 * in Hill, William G. 1975. TPB 8:117-126.
 */
double Hill_get_sigdsq(double y[], int twoNsmp) {
    double      sigdsq;

    if(twoNsmp == 0) {
        /* no bias correction */
        sigdsq = y[2] / (2.0 * y[0]);
    } else {
        /*
         * Sampling is like one generation of evolution,
         * with twoN = twoNsmp, and u = c = 0.
         */
        double      y2[DIM];

        onestep(y, y2, twoNsmp, 0.0, 0.0);
        sigdsq = y2[2] / (2.0 * y2[0]);
    }
    return sigdsq;
}

double Hill_exactLD(double c, double u, PopHist *ph, int twoNsmp, void *data) {
    HillData *hd = (HillData *) data;
    Hill_evolveDiscrete(hd->y, ph, c, u);
    return Hill_get_sigdsq(hd->y, twoNsmp);
}

/**
 * Evolve through all Epochs of population history, using step size h
 * and beginning with initial equilibrium. Calculation iterates Hill's
 * difference equation. The initial value of y is set using Hill's
 * equation for equilibrium.
 *
 * @param[out] y On return, y contains Hill's vector of moments.
 * @param[in] ph Describes the population's history. If ph contains
 * just one Epoch, y gets its equilibrium value.
 * @param[in] c Recombination rate; overrides values in ph.
 * @param[in] u Mutation rate.
 * @returns Always returns 0.
 */
int Hill_evolveDiscrete(double y[], PopHist * ph, double c, double u) {
    int         i;
    unsigned    nepoch = PopHist_nepoch(ph);

    /* set initial y to equilibrium */
    Hill_geteq(y, PopHist_twoN(ph, nepoch - 1), c, u);

#ifndef NDEBUG
    for(i = 0; i < DIM; ++i) {
        if(y[i] < 0.0) {
            printf("Hill_evolveDiscrete: negative y[%d]: %lg\n", i, y[i]);
        }
    }
#endif

#ifdef DEBUG
    int         status;

    status = pthread_mutex_lock(&outputLock);
    if(status)
        ERR(status, "lock output");

    assertFiniteArray(y, DIM, __FILE__, __LINE__);

    status = pthread_mutex_unlock(&outputLock);
    if(status)
        ERR(status, "unlock output");
#endif

    for(i = nepoch - 2; i >= 0; --i) {
        iterate(y, (int) PopHist_duration(ph, i), PopHist_twoN(ph, i), c, u);
    }

    return 0;
}

#ifndef NDEBUG
/*
 * Uses definitions on Hill's p 119.
 */
int getDRM_slow(double DRM[][DIM], double twoN, const double c,
				double u) {
    int         i;
    double      N = 0.5 * twoN; /* number of diploids individuals */
    double      m, m4, one_minus_c;
    double      A[DIM][DIM];
    double      RM[DIM];

    static const double P[DIM][DIM] = {
        {1.0, -0.5, 0.0, 0.0},
        {0.0, 2.5, -1.0, 0.0},
        {-1.0, -1.0, 1.5, 0.0},
        {0.0, 0.0, 0.0, 0.5}
    };

    /* Construct RM vector (recombination and mutation) */
    m = (1.0 - u);
    assert(m > 0.0);
    m4 = m * m * m * m;
    one_minus_c = 1.0 - c;
    assert(one_minus_c > 0.0);
    RM[0] = m4;
    RM[1] = one_minus_c * m4;
    RM[2] = one_minus_c * one_minus_c * m4;
    RM[3] = m * m;

    /*
     * Construct drift matrix, A = (I - P/N)*Diag(RM)
     * Set DRM equal to result.
     */
    memcpy(A, P, (size_t) (DIM * DIM * sizeof(A[0][0])));

    for(i = 0; i < DIM; ++i) {
        /* unroll inner loops for speed */
        A[i][0] /= -N;
        A[i][1] /= -N;
        A[i][2] /= -N;
        A[i][3] /= -N;
        A[i][i] += 1.0;

        /* matrix product A*Diag(RM) */
        DRM[i][0] = A[i][0] * RM[0];
        DRM[i][1] = A[i][1] * RM[1];
        DRM[i][2] = A[i][2] * RM[2];
        DRM[i][3] = A[i][3] * RM[3];
    }

    return 0;
}
#endif

#ifndef NDEBUG
/**
 * This version does not have unrolled loops.
 * Move the state vector forward one generation, using Hill's
 * recurrence (Eqns 5 and 8).
 *
 * @param[input] x is the initial state vector.
 * @param[output] y is the new state vector.
 */
int onestep_slow(const double x[], double y[], double twoN, double c,
				 double u) {
    int         i, j;

    /*double U = twoN*u; */
    double      DRM[DIM][DIM];

#ifdef DEBUG
    assertFiniteArray(x, DIM, __FILE__, __LINE__);
#endif

    /* DRM is matrix of effects of drift, recombination and mutation */
    getDRM_slow(DRM, twoN, c, u);

#ifdef DEBUG
    if(!matIsFinite(DIM, DRM)) {
        printf("onestep got nonfinite DRM matrix.\n");
        printf("  on entry: [%lg, %lg, %lg, %lg]\n", x[0], x[1], x[2], x[3]);
        printsqrmat("nonfinite DRM", DIM, DRM);
        exit(1);
    }
#endif

    /* begin with increment due to mutation (Hill's Eqn 8). */
    y[0] = 4.0 * u * x[3];
    y[1] = y[2] = 0.0;
    y[3] = 2.0 * u;

    /* Add the matrix-vector product, DRM*y (Hill's eqn 5) */
    for(i = 0; i < DIM; ++i) {
        for(j = 0; j < DIM; ++j)
            y[i] += DRM[i][j] * x[j];
    }

#ifdef DEBUG
    assertFiniteArray(y, DIM, __FILE__, __LINE__);
#endif

    return 0;
}
#endif
