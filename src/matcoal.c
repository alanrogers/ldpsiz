/**
 * @file matcoal.c
 * @brief Implement class MatCoal
 * @internal
 * Copyright (C) 2014 Alan R. Rogers
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * This file is released under the Internet Systems Consortium
 * License, which can be found in file "LICENSE".
 *
 * Alan R. Rogers, Department of Anthropology, University of Utah,
 * Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
 **/

#include "matcoal.h"
#include "misc.h"
#include "pophist.h"
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

/**
 * Class MatCoal: Matrix Coalescent (Wooding and Rogers 2001)
 *
 * The Matrix Coalescent is a re-formulation of the coalescent, a
 * stochastic process that is widely used within population
 * genetics.  The coalescent is a model of the genealogy of the
 * genes in a sample.  As one traces the genealogy backwards in
 * time, the number of lineages within it decreases by one each
 * time a branch point (or coalescent event) is encountered.
 * Beginning with a sample of n genes, the first coalescent event
 * leaves us with a sample of n-1, the second with a sample of
 * n-2, and so on until only a single gene is left, the common
 * ancestor of all the genes in the sample.  At any given time,
 * the process may be in any of n states.  A process in state i at
 * time t, where i is an element of (1,2,...,n), is one with i
 * distinct lineages in its genealogy at time t.  Here, time is
 * measured in generations backwards from the present.
 *
 * The matrix coalescent describes the same process using the
 * machinery of Markov chains.  The state of the process at time t
 * is described by a probability distribution: a vector of n
 * non-negative entries that sum to 1.  The i'th entry in the
 * vector (where i=1,2,...,n) is the probability that the
 * coalescent process is in state i at time t.  Although time is
 * measured in generations, it is treated as a continuous
 * variable.
 *
 * The virtue of this reformulation of the coalescent process is
 * that it deals gracefully with populations that vary in size.
 * For simplicity, the MatCoal class assumes that the history
 * population size can be divided into a series of epochs within
 * each of which the size is invariant.  The number of epochs,
 * however, is unlimited.
 *
 * The matrix coalescent does not deal gracefully with population
 * subdivision.  Or if it does, we have never figured out how.
 * Consequently, the MatCoal class assumes that the population
 * consists of a single randomly-mating deme.
 *
 * @author Alan Rogers
 * @date 18 November 2014
 **/

/**
 * Project vector forward in time assuming constant population
 * size.
 *
 * @param x initial vector.  On entry, x[i] should equal the
 * probability that there are i+1 separate lineages at time 0.
 * On return, x[i] equals the probability that there are i+1
 * lineages at time t.
 *
 * @param nSamples length of vector, which is equal to the sample
 * size.  Must equal the value of n that was used in constructing the
 * MatCoal object.
 *
 * @param v Number of time units to project the vector
 * forward.  Time is measured in units of twoN generations, where
 * twoN is the (constant) size of the population.  In other
 * words, v = g/twoN, where g is time in generations.
 *
 * This method projects a probability distribution backwards in
 * time using a continuous-time Markov Chain with transition rate
 * matrix A.  Matrix A has diagonal entries A[i][i] = -MatCoal_beta(i) and
 * super-diagonal entries A[i][i+1] = MatCoal_beta(i+1).  All the other
 * entries of A are zero.
 *
 * To project the matrix backwards in time, this method uses the
 * "Uniformization" algorithm that is described in chapter 8 of
 * Introduction to Numerical Solution of Markov Chains, by
 * W. J. Stewart.  In the documentation below, I refer to this book
 * as INSMC.
 *
 * The uniformization algorithm is an approximation to exp(A*v)*x,
 * for matrix A, vector x, and scalar v:
 *
 * Set g = the absolute value of the largest diagonal entry of A,
 * gv=g*v, P = I + A/g, and z = x*exp(-gv). Then approximately
 * @verbatim
 *
 * exp(A*v)*x ~= sum_{i=0}^K (gv)^i/i! * P^i * z
 *	    
 * = z + gv*P*z + gv*gv*P*P*z/2 + gv*gv*gv*P*P*P*z/(2*3) + etc
 * @endverbatim
 *
 * It is easy to calculate the optimal number of terms in advance,
 * making it possible to optimize this calculation by organizing it
 * in "Horner form".  For example, if 4 terms are retained, the
 * Horner version is
 * @verbatim
 *
 * z + gv*P(z + (gv*P/2)*(z + (gv*P/3)*z))
 * @endverbatim 
 *
 * If v is large, it is useful to carry out the calculation in
 * several steps.   The goal is to calculate x[v] = exp(A*v)*x[0].
 * This can be done in J steps as follows:
 * @verbatim
 *
 * h = v/J
 * x[  h] = exp(A*h)*x[0]
 * x[2*h] = exp(A*h)*x[h]
 * ...
 * x[J*h] = exp(A*h)*x[(J-1)*h]
 * @endverbatim
 * If the error in each step is eps then the total error is no
 * greater than J*eps.  Thus, the number of terms in the
 * approximation at each step should be chosen to keep the per-step
 * error within errTol/J, where errTol is some pre-determined
 * error tolerance.
 */

void MatCoal_project(unsigned nSamples, double *x, double v,
                     double betavec[nSamples], double errTol) {

    DPRINTF(("%s:%d %lu entry\n", __func__, __LINE__,
             (long unsigned) pthread_self()));

    if(v < 0.0)
        printf("%s:%d: v=%lg\n", __FILE__, __LINE__, v);
    assert(v >= 0.0);

	if(v == 0.0) {
        DPRINTF(("%s:%d %lu early return\n", __func__, __LINE__,
                 (long unsigned) pthread_self()));
		return;
    }

	assert(v > 0);

    // maxGV is the largest value such that exp(-maxGV) is doesn't underflow
    double maxGV = nextafter(-log(DBL_MIN), 0.0);    
	double g = betavec[nSamples-1]; /* abs(largest B[i][i]) */
    double gv = g*v;
    double y[nSamples];

    DPRINTF(("%s:%d %lu gv=%lf\n", __func__, __LINE__,
             (long unsigned) pthread_self(), gv));
	assert(g > 0);

	// If v is large, then v is broken into several steps.  This
	// paragraph determines the number, J, of steps.
	unsigned J = 1;
	if(gv > maxGV) {
		J = 1 + (unsigned) (gv/maxGV);
        errTol /= J;
        v /= J;
        gv = g*v;
    }

    DPRINTF(("%s:%d %lu doing %u steps\n", __func__, __LINE__,
             (long unsigned) pthread_self(), J));

	// Find K, which determines the number of terms in the
	// polynomial approximation.  This polynomial is of order K
	// and contains K+1 terms.  See INSMC, eqn 8.6, p 410.
	double s=1, p=1;
	double expMinusGV, rhs;

    errno = 0;
	expMinusGV = exp(-gv);
    if(errno) 
        eprintf("%s:%d: exp failed: %s\n", __FILE__, __LINE__,
                strerror(errno));
    assert(isnormal(expMinusGV));

	rhs = (1.0 - errTol)/expMinusGV;
	unsigned K=0;
	while(s < rhs){
		++K;
		p *= gv/K;
		s += p;
	}
    DPRINTF(("%s:%d %lu doing %u polynomial terms; J*K=%u\n",
             __func__, __LINE__,
             (long unsigned) pthread_self(), K, J*K));

	// Each pass through the loop below performs one step
	unsigned step;
	for(step=0; step < J; ++step) {
		unsigned i;

		// Initialize by setting y = x*exp(-g*v).  The initial
		// entry of x plays no role in the calculation and is
		// therefore ignored.
		for(i=1; i < nSamples; ++i)
			x[i] *= expMinusGV;

		// Copy x into y
		memcpy(y+1, x+1, (nSamples - 1) * sizeof(x[0]));

		// The next loop is the uniformization approximation
		// (INSMC eqn 8.4, p 410).  Each pass calculates y = x +
		// gv*P*y/i, where i = K..1.
		for(i=K; i > 0; --i) {
			double p0, p1;
			unsigned j;
			for(j=1; j < nSamples - 1; ++j) {
				p0 = 1.0 - betavec[j]/g;
				p1 = betavec[j+1]/g;
				y[j] = x[j] + (p0* y[j] + p1* y[j+1])*gv/i;
			}
			p0 = 1.0 - betavec[nSamples - 1]/g;
			y[nSamples - 1] = x[nSamples - 1] + p0* y[nSamples - 1]*gv/i;
		}

		// Copy y into x and calculate x[0]
		x[0] = 1;
		for(i=1; i < nSamples; ++i) {
			x[0] -= y[i];  // x[0] = 1 - sum of other elements
			x[i] = y[i];
		}
	}
    DPRINTF(("%s %lu exit\n", __func__, (long unsigned) pthread_self()));
}

/**
 * Project the probability vector forward to times
 * tvec[0],...,tvec[k-1].
 *
 * @param x  an nTimes X nSamples matrix, into which the results are
 * placed.  On return, x[i][j] holds the probability that the
 * coalescent process contains j+1 distinct lineages at time
 * tvec[i].
 *
 * @param tvec A k-vector whose entries give the times for which
 * results are to be calculated.
 **/
void MatCoal_project_multi(unsigned nTimes,  unsigned nSamples,
						   double x[nTimes][nSamples], double tvec[nTimes],
						   PopHist *ph, double errTol) {
    DPRINTF(("%s %lu entry\n", __func__, (long unsigned) pthread_self()));
    assert(nSamples != 0);

	// create initial vector.
	double w[nSamples];
	memset(w, 0, (nSamples-1)*sizeof(double));
	w[nSamples - 1] = 1.0;

    unsigned i;

    // betavec[i] is i(i+1)/2
    double betavec[nSamples];
    for(i=0; i<nSamples; ++i)
        betavec[i] = MatCoal_beta(i);

    int iepoch = 0;

	double t = 0;            // current time   
	double r = PopHist_duration(ph, iepoch); // time remaining in epoch
	unsigned tndx=0;         // index into tvec

	// Each pass through loop processes either the interval
	// between t and the next value of tvec or the interval
	// between t and the end of the current history epoch.
	while(tndx < nTimes) {

		if(!isfinite(PopHist_duration(ph, iepoch)) || tvec[tndx] <= t + r) {
			// Project to next tvec value (w/i current epoch)
			double step = tvec[tndx] - t;
			MatCoal_project(nSamples, w, step*PopHist_twoNinv(ph, iepoch),
                            betavec, errTol);

			for(i=0; i < nSamples; ++i)
                x[tndx][i] = w[i];
			r -= step;
			t = tvec[tndx];
			++tndx;
		}else {
			// Finish current epoch
			MatCoal_project(nSamples, w, r*PopHist_twoNinv(ph, iepoch),
                            betavec, errTol);
			t += r;
			++iepoch;
			r = PopHist_duration(ph, iepoch);
		}
	}
    DPRINTF(("%s %lu exit\n", __func__, (long unsigned) pthread_self()));
}

/**
 * Find the expected length of each coalescent interval within a
 * single epoch of population history, during which population size
 * was constant.
 *
 * This is done by integrating exp(A*t)*x from t=0 to
 * dt (which may be infinite).  Here, A is the transition rate matrix,
 * and x = * (0,..,0,1)' is the initial probability vector.
 *
 * @param nSamples The number of lineages sampled. Also the length of
 * vectors m, p0, and p1.
 *
 * @param p0 An array of doubles, representing the initial probability
 * vector at the beginning (tipward) end of the epoch. Value is not
 * changed.
 *
 * @param dt The length of the epoch in generations. Should equal
 * HUGE_VAL for the earliest history epoch.
 *
 * @param twoN Haploid population size during epoch.
 *
 * @param errTol Error tolerance.
 *
 * @param p1 An array of doubles of length nSamples. If dt is finite,
 * p1 will be set to the probability vector at the rootward end of the
 * epoch. If dt is infinite, p1 is not accessed.
 *
 * @param m On return m[i] holds the expected length (during the
 * epoch) of the coalescent interval with i+1 distinct lineages.  m[0]
 * is set equal to HUGE_VAL. 
 */
void MatCoal_integrate_epoch(unsigned nSamples,
                             double p0[nSamples],
                             double dt, double twoN, double errTol,
                             double p1[nSamples], double m[nSamples],
                             double betavec[nSamples]) {
    DPRINTF(("%s %lu entry\n", __func__, (long unsigned) pthread_self()));
    unsigned i;

    if(isfinite(dt)) {
        // Project probability across current epoch. Put result
        // in p1.
        memcpy(p1, p0, nSamples*sizeof(p0[0]));
        MatCoal_project(nSamples, p1, dt/twoN, betavec, errTol);

        // Set m equal to difference between probability vectors
        // at either end of current epoch.
        for(i=1; i<nSamples; ++i)
            m[i] = p1[i] - p0[i];

        // Form m = A^{-1} m
        for(i= nSamples-2; i > 0; --i) {
            m[i] += m[i+1];
            m[i+1] *= (-twoN/betavec[i+1]);
        }
        m[1] *= (-twoN/betavec[1]);
    }else{
        // Deal with infinite epoch at end of history.
        // The integral we need is -A^{-1} * x(nEpoch-1).

        // Form m = - A^{-1} x.
        m[nSamples-1] = p0[nSamples-1];
        for(i = nSamples-2; i > 0; --i) {
            m[i] = p0[i] + m[i+1];
            m[i+1] *= (twoN/betavec[i+1]);
        }
        m[1] *= (twoN / betavec[1]);
    }
    DPRINTF(("%s %lu exit\n", __func__, (long unsigned) pthread_self()));
}

/**
 * Find the expected length of each coalescent interval.
 *
 * This is done by integrating exp(A*t)*x from t=0 to
 * t=infinity.  Here, A is the transition rate matrix, and x =
 * (0,..,0,1)' is the initial probability vector.
 *
 * @param nSamples The number of lineages sampled. Also the length of
 * m and of each row of pr.
 *
 * @param m An nSamples-vector into which the answers will be
 * placed.  On return m[i] holds the expected length of
 * the coalescent interval during which there were i+1
 * distinct lineages.  m[0] is set equal to HUGE_VAL.
 *
 * The goal of this method is to calculate
 * @verbatim
 * F(0,infinity)*x(0)
 * @endverbatim
 * where F(a,b) is the definite integral of exp(A*z/N(t)) dz over
 * the interval (a,b), and 
 * @verbatim
 * x(0) = (0,...,0,1)'
 * @endverbatim
 * is the initial probability vector.
 *
 * It is useful to re-express this problem in terms of the
 * variable v(t) = t/N(t).  With this transformation, we have the
 * definite integral of N(z)*exp(A*v) dv over the interval
 * (b/N(b), a/N(a)).
 *
 * If the population's history is divided into P+1 epochs with
 * population sizes N0, N1, ..., NP, this integral can be
 * represented as a sum of smaller integrals.
 * @verbatim
 * G(0,infinity)*x(0) = G(0,v1)*x(0)
 *                      + G(v1,v2)*x(0) + ... + G(vP,infinity)*x(0)
 * @endverbatim
 * Here, interval (0,v1) corresponds to history epoch 0, (v1,v2) to
 * epoch 1, and (vP,infinity) to epoch P.  This decomposition of the
 * problem helps because population size within each epoch is
 * constant, so each of the integrals G(vi,vj) has a simple form.  For
 * epochs of finite length, 
 * @verbatim
 *                  vj
 *                  /
 * G(vi, vj)*x(0) = | Ni * exp(A*v) dv
 *                  /
 *                 vi
 *                = Ni*A^{-1} * (exp(A*vj) - exp(A*vi)) * x(0)
 *                = Ni*A^{-1} * ( x(vj) - x(vi) )
 * @endverbatim
 * For the final epoch, which has infinite length, this becomes
 * @verbatim
 *                  inf
 *                   /
 * G(vi, inf)*x(0) = | Ni * exp(A*v) dv
 *                   /
 *                   vi
 *                   = -Ni*A^{-1} * exp(A*vi)*x(0)
 * @endverbatim
 *
 * The algorithm used here first calls "project" to calculate the
 * probability vectors x(vi), (i=0,1,...), then subtracts pairs of
 * vectors, and finally applies Ni*A^{-1}.  This last step is easy
 * because 
 * @verbatim
 *                           [u2 + u3 + u4  u3 + u4   u4 ]'
 * A^{-1} * [u2, u3, u4]' =  [------------, -------, ----]
 *                           [     -a2         -a3    -a4]
 * @endverbatim
 * where -a2, -a3, and -a4 are the diagonal entries of A.
 *
 * @param m A pointer to a vector of n numbers of type double.
 * On return, m[i] = the expected length of the coalescent
 * interval that contains i+1 distinct lineages.  m[0]=infinity
 * because the coalescent interval containing a single lineage
 * lasts forever.
 *
 * @param nSamples The length of the vector.
 */
void MatCoal_integrate(unsigned nSamples, double m[nSamples], PopHist *ph,
                       double errTol) {
    DPRINTF(("%s %lu entry\n", __func__, (long unsigned) pthread_self()));
	unsigned i;
	unsigned nEpochs = PopHist_nepoch(ph); // number of epochs
	unsigned currEpoch;
	double twoN, dt;

	// pr[i][j]=Pr[j+1 lines of descent at tipward end of epoch i]
	double pr[nEpochs][nSamples]; 

	// Initial value, for tipward end of epoch zero, says that
	// there are n lines of descent with probability 1.
	memset(pr[0], 0, (nSamples-1)*sizeof(double));
	pr[0][nSamples-1] = 1.0;

	double y[nSamples];

    double betavec[nSamples];
    for(i=0; i<nSamples; ++i)
        betavec[i] = MatCoal_beta(i);

	// initialize m
	memset(m, 0, nSamples * sizeof(double));
	m[0] = HUGE_VAL;

	// Fill rows of pr matrix with probability vectors
	for(currEpoch=0; currEpoch < nEpochs; ++currEpoch)	{
		twoN = PopHist_twoN(ph, currEpoch);
		dt = PopHist_duration(ph, currEpoch);

        // Integrate current epoch, putting integral into y,
        // and setting pr[currEpoch+1] equal to probability vector
        // at end of current epoch.
        MatCoal_integrate_epoch(nSamples, pr[currEpoch],
                                dt, twoN, errTol,
                                pr[currEpoch+1], y,
                                betavec);

		// Add y to m
		for(i=1; i < nSamples; ++i)
			m[i] += y[i];
		
	}
    DPRINTF(("%s %lu exit\n", __func__, (long unsigned) pthread_self()));
}
