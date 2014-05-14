/**
 * @file strobeck.c
 * @author Alan R. Rogers
 * @brief Calculate sigma_d^2 using recurrence equations of Strobeck
 * and Morgan (1978), and Hudson (1985). 
 *
 * Strobeck and Morgan. 1978. The effect of intragenic recombination
 * on the number of alleles in a finite population. Genetics 88:
 * 829-844.
 *
 * It also uses formula A6 (p 630) from
 *
 * Hudson, Richard R. 1985.  The sampling distribution of linkage
 * disequilibrium under an infinite allele model without
 * selection. Genetics, 109(3): 611-631.
 *
 * Hudson showed that E[D^2] = Phi_AB - 2 Gamma_AB + Delta_AB
 *
 * The parameter sigma_d^2 was defined by
 *
 * Ohta, T. and Kimura, Motoo. 1971. Linkage disequilibrium between
 * two segregating nucleotide sites under the steady flux of mutations
 * in a finite population. Genetics 68(4): 571--580.
 *
 * They defined it like this
 *
 *   sigma_d^2 = E[D^2]/E[(1-J_A)(1-J_B)]
 *
 * where J_A = sum a_i^2, and J_B = sum b_i^2, where a_i and b_i are
 * allele frequencies of the current generation at loci A and B. The
 * denominator is thus the probability that two random genes from
 * locus A differ, and two genes drawn independently from locus B also
 * differ. 
 *
 * It can also be shown (I haven't found this in the literature) that
 *
 * E[(1-J_A)(1-J_B)] = 1 - PhiA - PhiB + DeltaAB
 *
 * The code below uses the formulas of Strobeck and Morgan to
 * calculate the 5 identity coefficients, then uses these to calculate
 * sigma_d^2. 
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include "pophist.h"
#include "misc.h"
#include "model.h"
#include "string.h"
#include "strobeck.h"

#define DIM 5

static const char *stateLbl[DIM] = { "phiA",
    "phiB",
    "phiAB",
    "gammaAB",
    "deltaAB"
};

typedef struct {
    double      twoN;           /* haploid pop size */
    double      c;              /* recombination rate per generation */
    double      u;              /* mutation rate per generation */
} param_struct;

typedef struct {
    unsigned    ydim;
    double      y[DIM];
} StrobeckData;

double Strobeck_exactLD(double c, double u, PopHist *ph, int twoNsmp,
                        void *data);

static int  onestep(const double y1[], double y2[], double twoN, double c,
                    double u);
static int  iterate(double y[], int t, double twoN, double c, double u);
static void geteq_c_mid(double x[DIM], double twoN, double c,
                        double u, double fourNu);

void       *StrobeckData_new(void) {
    StrobeckData *hd = malloc(sizeof(StrobeckData));
    checkmem(hd, __FILE__, __LINE__);
    memset(hd, 0, sizeof(StrobeckData));
    hd->ydim = DIM;

    return hd;
}

/** i'th state value */
double StrobeckData_stateVal(void *vdata, unsigned i) {
    StrobeckData *data = vdata;

    return data->y[i];
}

void StrobeckData_free(void *p) {
    StrobeckData *data = p;

    free(data);
}

Model      *Model_allocStrobeck(int twoNsmp) {
    Model      *model = malloc(sizeof(Model));

    checkmem(model, __FILE__, __LINE__);

    snprintf(model->lbl, sizeof(model->lbl), "Strobeck");
    model->twoNsmp = twoNsmp;

    model->ld = Strobeck_sigdsq;
    model->exactLD = Strobeck_exactLD;
    model->ldEq = Strobeck_sigdsqEq;
    model->stateDim = Strobeck_stateDim;
    model->stateLbl = Strobeck_stateLbl;
    model->stateVal = StrobeckData_stateVal;
    model->newState = StrobeckData_new;
    model->freeState = StrobeckData_free;

    return model;
};

/**
 * Set state vector to equilibrium for epoch "whichEpoch" and then use it
 * to calculate sigmdsq.
 */
double Strobeck_sigdsqEq(double c, double u, PopHist * ph,
                         unsigned whichEpoch, int twoNsmp, void *vdata) {

    StrobeckData *data = (StrobeckData *) vdata;

    /* y is equilibrium for appropriate epoch */
    Strobeck_geteq(data->y, PopHist_twoN(ph, whichEpoch), c, u);

    return Strobeck_get_sigdsq(data->y, twoNsmp);
}

/**
 * Set state vector to initial equilibrium, then evolve it through the
 * population history, and then use this state vector to calculate sigmdsq.
 */
double Strobeck_sigdsq(ODE * ode, double c, double u, PopHist * ph,
                       int twoNsmp) {
    StrobeckData *data = ODE_state(ode);
    unsigned    nepoch = PopHist_nepoch(ph);

    /* Initial y is equilibrium for earliest epoch */
    Strobeck_geteq(data->y, PopHist_twoN(ph, nepoch - 1), c, u);

    /* evolve through population history */
    ODE_evolve(ode, data->y, data->ydim, c, u, ph, Strobeck_dydt, !VERBOSE);

    return Strobeck_get_sigdsq(data->y, twoNsmp);
}

/* dimension of state vector */
size_t Strobeck_stateDim(void) {
    return DIM;
}

const char *Strobeck_stateLbl(unsigned i) {
    return stateLbl[i];
}

/**
 *
 * This function calculates the difference, from one generation to the
 * next, for each of the 5 identity coefficients. It is based on
 * equation 3 of Strobeck and Morgan. 
 *
 * y[0] = phiA
 * y[1] = phiB
 * y[2] = phiAB
 * y[3] = gammaAB
 * y[4] = deltaAB
 *
 * @param[in] y1 is the current vector of identity coefficients.
 * @param[out] f will contain the derivatives of y[], as approximated by
 * the difference equations of Strobeck and Morgan.
 * @param[in] params points to a structure of type param_struct.
 */
int
Strobeck_dydt(double t_notused, const double yy1[], double f[], void *params) 
{
    int         i;
    double      y2[DIM];
    struct dydt_params *par = (struct dydt_params *) params;

#ifdef DEBUG
    assertFiniteArray(yy1, DIM, __FILE__, __LINE__);
#endif

    onestep(yy1, y2, par->twoN, par->c, par->u);

#ifdef DEBUG
    assertFiniteArray(y2, DIM, __FILE__, __LINE__);
#endif

    /* f is the difference between new and old state vectors. */
    for(i = 0; i < DIM; ++i)
        f[i] = y2[i] - yy1[i];

#ifdef DEBUG
    assertFiniteArray(f, DIM, __FILE__, __LINE__);
#endif

    return GSL_SUCCESS;
}

#if 0
static void geteq_c_small(double x[DIM], double fourNu);
static void geteq_c_large(double x[DIM], double fourNu);
void geteq_c_small(double x[DIM], double fourNu) {
    x[0] = x[1] = 1.0 / (1.0 + fourNu);
    x[2] = 1.0 / (1.0 + 2.0 * fourNu);
    x[3] = 3.0 + 5 * fourNu;    /* numerator of GammaAB */
    x[3] /= (1.0 + fourNu) * (1 + 2 * fourNu) * (3.0 + 2 * fourNu);
    x[4] = 9.0 + 18 * fourNu + 4 * fourNu * fourNu; /* numerator DeltaAB */
    x[4] /= (1.0 + fourNu) * (1 + 2 * fourNu) * (3.0 + fourNu) * (3.0 +
                                                                  2 * fourNu);
}

/*
 * Strobeck and Morgan (p 834) give this as the equilibrium when
 * linkage is loose (c >> u). It always implies, however, that Dsq=0,
 * using Hudson's formula, Dsq = x[2] - 2.0*x[3] + x[4], which is
 * always zero. This makes no sense. I'm using the equilibrium formula
 * for intermediate c even when c is large. This present function is
 * not used.
 */
void geteq_c_large(double x[DIM], double fourNu) {
    x[0] = x[1] = 1.0 / (1.0 + fourNu);
    x[2] = x[3] = x[4] = x[0] * x[0];
}
#endif

void geteq_c_mid(double x[DIM], double twoN, double c, double u,
                 double fourNu) {
    double      den;
    double      uu = u * u, cc = c * c;

    x[0] = x[1] = 1.0 / (1.0 + fourNu);

    /*************************************************
    phiAB, gammaAB, and deltaAB have the same denominator. den
    expresses it in Horner form. I used Maxima's "horner" command
    to do the conversion. Here it is before conversion to Horner
    form:
    
    den: (2*u*twoN+1)*(32*u^3*twoN^3+24*c*u^2*twoN^3+4*c^2*u*twoN^3
                             +80*u^2*twoN^2+38*c*u*twoN^2+2*c^2*twoN^2
                             +54*u*twoN+13*c*twoN+9);
    *******************************************************/
    den = twoN * (twoN * (twoN * (uu * (u * (64 * u + 48 * c) + 8 * cc) * twoN
                                  + u * (u * (192 * u + 100 * c) + 8 * cc))
                          + u * (188 * u + 64 * c) + 2 * cc)
                  + 72 * u + 13 * c) + 9;

    /*************************************************
    PhiAB. Numerator and denominator in Horner form. Here is the
    numerator before conversion to Horner:  
    
    num: 16*u^3*twoN^3+4*c*u^2*twoN^3+44*u^2*twoN^2+12*c*u*twoN^2+2*c^2*twoN^2
             +36*u*twoN+13*c*twoN+9;
    *******************************************************/
    x[2] =
        twoN * (twoN *
                (uu * (16 * u + 4 * c) * twoN + u * (44 * u + 12 * c) +
                 2 * cc) + 36 * u + 13 * c) + 9;
    x[2] /= den;

    /******************************************************
    gammaAB: numerator and denominator in Horner form
    Numerator before conversion to Horner:  

    num: 20*u^2*twoN^2 + 12*c*u*twoN^2 + 2*c^2*twoN^2 + 36*u*twoN
         + 13*c*twoN + 9;
    *******************************************************/
    x[3] =
        twoN * ((u * (20 * u + 12 * c) + 2 * cc) * twoN + 36 * u + 13 * c) +
        9;
    x[3] /= den;

    /******************************************************
    deltaAB: numerator and denominator in Horner form. Numerator
    before conversion to horner:

    num: 16*u^2*twoN^2+12*c*u*twoN^2+2*c^2*twoN^2+36*u*twoN+13*c*twoN+9;
    *******************************************************/
    x[4] =
        twoN * ((u * (16 * u + 12 * c) + 2 * cc) * twoN + 36 * u + 13 * c) +
        9;
    x[4] /= den;
}

/**
 * Fill x vector with equilibrium values given on pp 833-834 of
 * Strobeck and Morgan.
 *
 * Strobeck and Morgan separate 3 cases, for c << u, c >> u, and c ~=
 * u. I did some experiments to figure out how to bound these ranges,
 * and decided never to use the formula for large c.
 */
int Strobeck_geteq(double x[], double twoN, double c, double u) {
    double      fourNu = 2.0 * twoN * u;

    /*
     * The formula for intermediate recombination always seems best,
     * so I'm ignoring the other two.
     */
    geteq_c_mid(x, twoN, c, u, fourNu);

    return 0;
}

/**
 * Calculate sigma_d^2 from state vector.
 *
 * The formulas involved are in the comment at the top of this file.
 *
 * @param[in] y State vector, whose entries represent PhiA, PhiB,
 * PhiAB, GammaAB, and DeltaAB.
 * @param[in] n Number of gene copies in sample. If n==0, then no
 * bias correction is done.
 */
double Strobeck_get_sigdsq(double y[], int twoNsmp) {
    double      sigdsq;

    if(twoNsmp == 0) {
        /* without bias correction */
        sigdsq = (y[2] - 2.0 * y[3] + y[4]) / (1.0 - y[0] - y[1] + y[4]);
    } else {

#if 0
        /* Hudson's (1985, p 631) bias correction */
        double      a = (twoNsmp - 1.0) / twoNsmp;
        double      b = (twoNsmp - 2.0) / twoNsmp,
            double c = (twoNsmp - 3.0) / twoNsmp;
        double      phiA = 1.0 / twoNsmp + a * y[0];
        double      phiB = 1.0 / twoNsmp + a * y[1];
        double      phiAB = 1.0 / twoNsmp + a * y[2];
        double      gammaAB = a * b * y[3] + 2.0 * a * y[0] / twoNsmp
            + a * y[2] / twoNsmp + 1.0 / (twoNsmp * twoNsmp);
        double      deltaAB =
            a * b * c * y[4] + 2 * a * b * (y[0] + 2.0 * y[3]) / twoNsmp +
            2 * a * (2 * y[0] + y[2]) / (twoNsmp * twoNsmp) +
            1.0 / (twoNsmp * twoNsmp);
        sigdsq =
            (phiAB - 2.0 * gammaAB + deltaAB) / (1.0 - phiA - phiB + deltaAB);
#else
        /*
         * Sampling is like one generation of evolution,
         * with twoN = twoNsmp, and u = c = 0.
         */
        double      y2[DIM];

        onestep(y, y2, twoNsmp, 0.0, 0.0);
        sigdsq =
            (y2[2] - 2.0 * y2[3] + y2[4]) / (1.0 - y2[0] - y2[1] + y2[4]);
#endif
    }

    return sigdsq;
}

/**
 * Move the state vector forward one generation, using recurrence
 * equations of Strobeck and Morgan (1978, see Eqn 3, p 833).
 *
 * @param[input] yy1 is the initial state vector.
 * @param[output] y2 is the new state vector.
 */
static int onestep(const double yy1[], double y2[], double twoN, double c,
                   double u) {

    double      twoNsq = twoN * twoN;
    double      twoNinv = 1.0 / twoN;
    double      omc = 1.0 - c;
    double      omu = 1.0 - u;
    double      omusq = omu * omu;
    double      omu4th = omusq * omusq;
    double      csq = c * c;

    /* gamma and delta are defined on p 831 of Strobeck and Morgan */
    double      gam = 1.0 + (twoN - 1.0) * (yy1[0] + yy1[1] + yy1[2])
        + (twoN - 1.0) * (twoN - 2.0) * yy1[3];

    gam /= twoNsq;

    double      delta = twoN + 2.0 * (twoN - 1.0) * (yy1[0] + yy1[1] + yy1[2])
        + (twoN - 1.0) * (twoN - 2.0) * (yy1[0] + yy1[1] + 4.0 * yy1[3])
        + (twoN - 1.0) * (twoN - 2.0) * (twoN - 3.0) * yy1[4];

    delta /= twoN * twoN * twoN;

    /*
     * Recurrences defined on p 833 of Strobeck and Morgan.
     */
    y2[0] = omusq * (twoNinv + (1.0 - twoNinv) * yy1[0]);
    y2[1] = omusq * (twoNinv + (1.0 - twoNinv) * yy1[1]);
    y2[2] = omu4th * (omc * omc * (twoNinv + (1.0 - twoNinv) * yy1[2])
                      + 2 * c * omc * gam + csq * delta);
    y2[3] = omu4th * (omc * gam + c * delta);
    y2[4] = omu4th * delta;

    return 0;
}

/**
 * Iterate the difference equation across t generations.
 */
static int iterate(double y[], int t, double twoN, double c, double u) {
    double      y2[DIM];

    while(t-- > 0) {
        onestep(y, y2, twoN, c, u);
        memcpy(y, y2, sizeof(y2));
    }
    return 0;
}

double Strobeck_exactLD(double c, double u, PopHist *ph, int twoNsmp,
                        void *data) {
    StrobeckData *sd = (StrobeckData *) data;

    Strobeck_evolveDiscrete(sd->y, ph, c, u);

    return Strobeck_get_sigdsq(sd->y, twoNsmp);
}

/**
 * Evolve through all epochs of population history, using step size h
 * and beginning with initial equilibrium. Calculation iterates 
 * difference equation. The initial value of y is set to equilibrium
 * equation for equilibrium.
 *
 * @param[out] y On return, y contains Hill's vector of moments.
 * @param[in] ph Describes the population's history. If ph contains
 * just one epoch, y gets its equilibrium value.
 * @param[in] c Recombination rate; overrides values in ph.
 * @param[in] u Mutation rate.
 * @returns Always returns 0.
 */
int Strobeck_evolveDiscrete(double y[], PopHist * ph, double c, double u) {
    int         i;
    unsigned    nepoch = PopHist_nepoch(ph);

    /* set initial y to equilibrium */
    Strobeck_geteq(y, PopHist_twoN(ph, nepoch - 1), c, u);

#ifdef DEBUG
    assertFiniteArray(y, DIM, __FILE__, __LINE__);
#endif

    for(i = PopHist_nepoch(ph) - 2; i >= 0; --i) {
        iterate(y, (int) PopHist_duration(ph, i), PopHist_twoN(ph, i), c, u);
    }

    return 0;
}

#if 0
int void    jac(double t, const double y[], double *dfdy,
                double dfdt[], void *params);

int void
jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    double      mu = *(double *) params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, -2.0 * mu * y[0] * y[1] - 1.0);
    gsl_matrix_set(m, 1, 1, -mu * (y[0] * y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}
#endif
