/**
 * @file em.c
 * @author Alan R. Rogers
 * @brief Estimate D from partially phased diploid data using the EM
 * algorithm.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "em.h"
#include "misc.h"

double      negLnL(double D, void *data);
double      negLnL_df(double D, void *data);
void        negLnL_fdf(double D, double *f, double *df, void *data);
void        negLnL_fddf(double D, double *f, double *df, double *ddf,
                        void *data);
void        negLnL_df12(double D, double *df, double *ddf, void *data);
static inline double get_w(double D, DsqData * dd);
static inline double getD(double w, DsqData * dd);
static inline double next_w(double w, DsqData * dd);
double      find_D(double D0, DsqData * dd);
static inline double initialD(DsqData * dd);

/** Print an object of type DsqData */
void DsqData_print(DsqData * dd, const char *file, int line, FILE * fp) {
    fprintf(fp, "%s:%d: DsqData is:\n", file, line);
    fprintf(fp, "  %-8s=%lg\n", "tol", dd->tol);
    fprintf(fp, "  %-8s=%lg\n", "alpha", dd->alpha);
    fprintf(fp, "  %-8s=%lg\n", "beta", dd->beta);
    fprintf(fp, "  %-8s=%lg\n", "loD", dd->loD);
    fprintf(fp, "  %-8s=%lg\n", "hiD", dd->hiD);
    fprintf(fp, "  %-8s=%lg\n", "px", dd->px);
    fprintf(fp, "  %-8s=%lg\n", "py", dd->py);
    fprintf(fp, "  %-8s=%u\n", "nGtype", dd->nGtype);
    fprintf(fp, "  %-8s=%u\n", "nUnphased", dd->nUnphased);
    fprintf(fp, "  %-8s=[%u, %u, %u, %u]\n", "nGam",
            dd->nGam[0], dd->nGam[1], dd->nGam[2], dd->nGam[3]);
}

/** Make DsqData object consistent with its current allele freqs.*/
void DsqData_reset(DsqData * dd) {
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1 - px;
    double      qy = 1 - py;

    dd->alpha = px * qx * py * qy;
    dd->beta = px + py - 2.0 * px * py;
    dd->loD = loD(px, py, dd->nGam);
    dd->hiD = hiD(px, py, dd->nGam);
}

/**
 * Return value of w implied by D, and the parameters in d.
 */
static inline double get_w(double D, DsqData * dd) {
    double      Z = dd->alpha - D * (dd->beta - D);
    double      w = (D + Z) / (D + 2.0 * Z);

    if(w < 0.0)
        return 0.0;
    if(w > 1.0)
        return 1.0;
    return w;
}

/**
 * Negative of log Likelihood of given value of D.
 */
double negLnL(double D, void *data) {
    DsqData    *dd = (DsqData *) data;

    double      x = 0.0, lnL = 0.0;
    double      w = get_w(D, dd);
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1.0 - dd->px;
    double      qy = 1.0 - dd->py;

    double      g[4] = { qx * qy + D,
        qx * py - D,
        px * qy - D,
        px * py + D
    };

    x = dd->nGam[0] + w * dd->nUnphased;
    if(x > 0)
        lnL += x * log(g[0]);

    x = dd->nGam[1] + (1 - w) * dd->nUnphased;
    if(x > 0)
        lnL += x * log(g[1]);

    x = dd->nGam[2] + (1 - w) * dd->nUnphased;
    if(x > 0)
        lnL += x * log(g[2]);

    x = dd->nGam[3] + w * dd->nUnphased;
    if(x > 0)
        lnL += x * log(g[3]);

    return -lnL;
}

/**
 * Negative of log Likelihood of given value of D. This version
 * calculates both the function value (placed in *f) and its
 * derivative (placed in *df).
 *
 */
void negLnL_fdf(double D, double *f, double *df, void *data) {
    DsqData    *dd = (DsqData *) data;

    double      x = 0.0;
    double      w = get_w(D, dd);
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1.0 - dd->px;
    double      qy = 1.0 - dd->py;

    double      g[4] = { qx * qy + D,
        qx * py - D,
        px * qy - D,
        px * py + D
    };

    *f = *df = 0.0;

    x = dd->nGam[0] + w * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[0]);
        *df += x / g[0];
    }

    x = dd->nGam[1] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[1]);
        *df -= x / g[1];
    }

    x = dd->nGam[2] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[2]);
        *df -= x / g[2];
    }

    x = dd->nGam[3] + w * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[3]);
        *df += x / g[3];
    }

    *f = -*f;
    *df = -*df;

    return;
}

/**
 * Negative of log Likelihood of given value of D. This version
 * returns only the 1st derivative.
 *
 */
double negLnL_df(double D, void *data) {
    DsqData    *dd = (DsqData *) data;

    double      x = 0.0;
    double      w = get_w(D, dd);
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1.0 - dd->px;
    double      qy = 1.0 - dd->py;

    double      g[4] = { qx * qy + D,
        qx * py - D,
        px * qy - D,
        px * py + D
    };

    double      df = 0.0;

    x = dd->nGam[0] + w * dd->nUnphased;
    if(x > 0) {
        df += x / g[0];
    }

    x = dd->nGam[1] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        df -= x / g[1];
    }

    x = dd->nGam[2] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        df -= x / g[2];
    }

    x = dd->nGam[3] + w * dd->nUnphased;
    if(x > 0) {
        df += x / g[3];
    }

    return -df;
}

/**
 * Negative of log Likelihood of given value of D. This version
 * calculates, the function value (placed in *f), its 1st
 * derivative (placed in *df), and its 2nd derivative (in ddf).
 *
 */
void negLnL_fddf(double D, double *f, double *df, double *ddf, void *data) {
    DsqData    *dd = (DsqData *) data;

    double      x = 0.0;
    double      w = get_w(D, dd);
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1.0 - dd->px;
    double      qy = 1.0 - dd->py;

    double      g[4] = { qx * qy + D,
        qx * py - D,
        px * qy - D,
        px * py + D
    };

    *f = *df = *ddf = 0.0;

    x = dd->nGam[0] + w * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[0]);
        *df += x / g[0];
        *ddf -= x / (g[0] * g[0]);
    }

    x = dd->nGam[1] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[1]);
        *df -= x / g[1];
        *ddf -= x / (g[1] * g[1]);
    }

    x = dd->nGam[2] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[2]);
        *df -= x / g[2];
        *ddf -= x / (g[2] * g[2]);
    }

    x = dd->nGam[3] + w * dd->nUnphased;
    if(x > 0) {
        *f += x * log(g[3]);
        *df += x / g[3];
        *ddf -= x / (g[3] * g[3]);
    }

    *f = -*f;
    *df = -*df;
    *ddf = -*ddf;

    return;
}

/**
 * Negative of log Likelihood of given value of D. This version
 * calculates the 1st derivative (placed in *df), and the 2nd
 * derivative (in ddf).
 *
 */
void negLnL_df12(double D, double *df, double *ddf, void *data) {
    DsqData    *dd = (DsqData *) data;

    double      x = 0.0;
    double      w = get_w(D, dd);
    double      px = dd->px;
    double      py = dd->py;
    double      qx = 1.0 - dd->px;
    double      qy = 1.0 - dd->py;

    double      g[4] = { qx * qy + D,
        qx * py - D,
        px * qy - D,
        px * py + D
    };

    *df = *ddf = 0.0;

    x = dd->nGam[0] + w * dd->nUnphased;
    if(x > 0) {
        *df += x / g[0];
        *ddf -= x / (g[0] * g[0]);
    }

    x = dd->nGam[1] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *df -= x / g[1];
        *ddf -= x / (g[1] * g[1]);
    }

    x = dd->nGam[2] + (1 - w) * dd->nUnphased;
    if(x > 0) {
        *df -= x / g[2];
        *ddf -= x / (g[2] * g[2]);
    }

    x = dd->nGam[3] + w * dd->nUnphased;
    if(x > 0) {
        *df += x / g[3];
        *ddf -= x / (g[3] * g[3]);
    }

    *df = -*df;
    *ddf = -*ddf;

    return;
}

/** Lowest feasible value of D, given allele frequencies */
double loD(double pA, double pB, unsigned *nGam) {
    double      qA = 1.0 - pA;
    double      qB = 1.0 - pB;
    double      lo = fmax(-pA * pB, -qA * qB);

    /* prevent log function from generating NaNs */
    if(nGam[0] || nGam[3])
        lo += DBL_EPSILON;

    return lo;
}

/** Highest feasible value of D, given allele frequencies */
double hiD(double pA, double pB, unsigned *nGam) {
    double      qA = 1.0 - pA;
    double      qB = 1.0 - pB;
    double      hi = fmin(pA * qB, qA * pB);

    /* prevent log function from generating NaNs */
    if(nGam[1] || nGam[2])
        hi -= DBL_EPSILON;

    return hi;
}

#define GO_RIGHT(F,A,X,B)   ((X) + (F)*((B)-(X)))
#define GO_LEFT(F,A,X,B) ((X) - (F)*((X)-(A)))

/**
 * This algorithm finds the minimum of function negLnL within the
 * range [a,b], where a = dd->loD, and b=dd->hiD.
 *
 * The algorithm uses Newton steps when these go downhill, and
 * otherwise bisects in the downhill direction.
 *
 * When Newton fails, the next step is in the downhill direction, as
 * indicated by the derivative. If the derivative at the boundary has
 * the same sign as that at the current point, the algorithm bets that
 * the function is monotonic and takes a big step. Otherwise, there
 * must be an intermediate minimum, so the algorithm bisects.
 *
 * @param[out] D Lewontin's original measure of linkage disequilibrium,
 *
 * @param[in] dd pointer to object of type DsqData
 *
 * @returns the integer 0
 */
int minimize1D(double *D, DsqData * dd) {

    double      x1, x0, df, ddf, dfa, dfb, diff;
    double      a = dd->loD, b = dd->hiD;
    unsigned    ok;             /* did Newton step work? */

    dfa = negLnL_df(a, dd);
    dfb = negLnL_df(b, dd);

    x0 = initialD(dd);          /* initial guess */

    do {
        assert(a <= x0);
        assert(x0 <= b);

        negLnL_df12(x0, &df, &ddf, dd);

        /* try Newton step */
        ok = 0;
        if(ddf > 0.0) {
            x1 = x0 - df / ddf;
            if(x1 >= a && x1 <= b)
                ok = 1;         /* worked */
        }

        if(!ok) {
            /*
             * Newton failed: use modified bisect.
             *
             * If df has same sign as derivative (dfa or dfb) at
             * boundary in the downhill direction, then bet that the
             * function is monotone and take a big step. Otherwise
             * there is an interior minimum, so bisect.
             */
            if(df < 0.0) {      /* downhill is to the right */
                a = x0;
                dfa = df;
                if(dfb < 0.0)
                    x1 = GO_RIGHT(0.8, a, x0, b);
                else
                    x1 = GO_RIGHT(0.5, a, x0, b);
            } else {            /* assume downhill is to the left */
                b = x0;
                dfb = df;
                if(dfa > 0.0)
                    x1 = GO_LEFT(0.8, a, x0, b);
                else
                    x1 = GO_LEFT(0.5, a, x0, b);
            }
        }

        diff = fabs(x1 - x0);
        x0 = x1;
    } while(diff > dd->tol);

    assert((fabs(x1 - dd->loD) >= dd->tol
            && fabs(x1 - dd->hiD) >= dd->tol && ddf > 0.0)
           || (fabs(x1 - dd->loD) <= dd->tol && df > 0.0)
           || (fabs(x1 - dd->hiD) <= dd->tol && df < 0.0));

    *D = x1;

    return 0;
}

/**
 * Return value of D implied by given value of w.
 *
 * @param[in] w Current guess as to fraction of unphased
 * double-heterozygotes with two-locus genotype 11/00 rather than
 * 10/01. 
 *
 * @param[in] dd pointer to object of type DsqData
 *
 * @returns improved estimate of w
 */
static inline double getD(double w, DsqData * dd) {
    double      sxy = dd->nGam[3] + (dd->nUnphased) * w;
    unsigned    twoN = 2 * (dd->nGtype);
    double      D = (sxy - twoN * dd->px * dd->py) / (twoN - 1);

    return D;
}

/**
 * Initialize DsqData object from current data, assuming that
 * half of unphased heterozygotes have 2-locus genotype 11/00 and the
 * other half have 10/01.
 */
static inline double initialD(DsqData * dd) {
    double      sxy = dd->nGam[3] + (dd->nUnphased) * 0.5;
    unsigned    twoN = 2 * (dd->nGtype);
    double      D = (sxy - twoN * dd->px * dd->py) / (twoN - 1);

    if(D < dd->loD)
        D = dd->loD + dd->tol;
    else if(D > dd->hiD)
        D = dd->hiD - dd->tol;

    return D;
}

/**
 * One step of EM algorithm.
 *
 * @param[in] w Current guess as to fraction of unphased
 * double-heterozygotes with two-locus genotype 11/00 rather than
 * 10/01. 
 *
 * @param[in] dd pointer to object of type DsqData
 *
 * @returns improved estimate of w
 */
static inline double next_w(double w, DsqData * dd) {
    double      D = getD(w, dd);

    return get_w(D, dd);
}

/**
 * Ad hoc iterations to find D. 
 * @param[in] D0 Starting value for EM algorithm.
 * @param[in] d  Points to a structure containing data.
 */
double find_D(double D0, DsqData * dd) {
    unsigned    itr = 0;
    const static unsigned maxIterations = 100;

    double      w, D = D0;

    /* iterate until D ~ D0 */
    do {
        ++itr;
        D0 = D;
        w = get_w(D0, dd);
        D = getD(w, dd);
    } while(fabs(D - D0) > dd->tol && ++itr < maxIterations);

    if(itr == maxIterations)
        fprintf(stderr, "Warning@%s:%d: find_D did not converge\n",
                __FILE__, __LINE__);

    return D;
}
