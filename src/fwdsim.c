/**
 * Forward simulation with two loci, mutation, recombination, and drift.
 */
/*
 * Internet Systems Consortium License
 * 
 * Copyright (c) 2014, Alan R. Rogers <rogers@anthro.utah.edu>
 * 
 * Permission to use, copy, modify, and/or distribute this software for
 * any purpose with or without fee is hereby granted, provided that the
 * above copyright notice and this permission notice appear in all
 * copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 * AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
 * TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "misc.h"

int main(int argc, char **argv) {

    long        i, rep, nreps, npoly, gen, maxgen;
    unsigned    y[4];
    unsigned    twoN;
    time_t      currtime = time(NULL);
    double      c, pA, pB, x[4], D, Dsq, pqpq, Dsqsum, pqpqsum;
    double      u, omu, omusq, uomu, usq;
    gsl_rng    *rng;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned) currtime);

    maxgen = 10000;
    nreps = 10000;
    twoN = 5000;
    c = 0.001;                  /* recombination rate */
    u = 1e-4;                   /* mutation rate */
    omu = 1.0 - u;              /* one minus u */
    omusq = omu * omu;
    uomu = u * omu;
    usq = u * u;

    npoly = 0;
    Dsqsum = pqpqsum = 0.0;

    printf("2N=%d c=%lf u=%lf\n", twoN, c, u);
    /*    printf("%5s %6s %6s %6s\n", "rep", "pA", "pB", "rsq"); */
    for(rep = 0; rep < nreps; ++rep) {
        x[0] = 1.0 / twoN;      /* freq AB */
        x[1] = 0.0;             /* freq Ab */
        x[2] = 0.5 - 1.0 / twoN;    /* freq aB */
        x[3] = 0.5;             /* freq ab */

        for(gen = 0; gen < maxgen; ++gen) {
            /* mutate */
            x[0] = x[0] * omusq + uomu * (x[1] + x[2]) + usq * x[3];
            x[1] = x[1] * omusq + uomu * (x[0] + x[3]) + usq * x[2];
            x[2] = x[2] * omusq + uomu * (x[0] + x[3]) + usq * x[1];
            x[3] = x[3] * omusq + uomu * (x[1] + x[2]) + usq * x[0];

            /* recombine */
            pA = x[0] + x[1];
            pB = x[0] + x[2];
            x[0] = x[0] * (1 - c) + c * pA * pB;
            x[1] = x[1] * (1 - c) + c * pA * (1 - pB);
            x[2] = x[2] * (1 - c) + c * (1 - pA) * pB;
            x[3] = x[3] * (1 - c) + c * (1 - pA) * (1 - pB);

            gsl_ran_multinomial(rng, 4u, twoN, x, y);
            for(i = 0; i < 4; ++i)
                x[i] = y[i] / (double) twoN;
            pA = x[0] + x[1];
            pB = x[0] + x[2];
            pqpq = pA * (1 - pA) * pB * (1 - pB);
            if(Dbl_near(pqpq, 0.0))
                break;
        }

        if(Dbl_near(pqpq, 0.0))
            continue;

        ++npoly;
        pqpqsum += pqpq;
        D = x[0] * x[3] - x[1] * x[2];
        Dsq = D * D;
        Dsqsum += Dsq;
        /* printf("%5ld %8.5lg %8.5lg\n", rep, Dsq, pqpq); */
    }
    printf("sigdsq(%ld) = %lf\n", npoly, Dsqsum / pqpqsum);

    gsl_rng_free(rng);
    return 0;
}
