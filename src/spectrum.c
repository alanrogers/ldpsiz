/**
 * @file spectrum.c
 * @author Alan R. Rogers
 * @brief Manipulate site frequency spectrum.
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "spectrum.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/// Site frequency spectrum
struct Spectrum {
	unsigned nSamples, nObs;
	unsigned long *v;
};

double Spectrum_value(Spectrum *s, unsigned i) {
    if(s->nObs == 0)
        eprintf("%s:%s:%d: nObs==0", __FILE__, __func__, __LINE__);
    return s->v[i]/((double) s->nObs);
}

Spectrum *Spectrum_new(unsigned nSamples) {
    Spectrum *spec = malloc(sizeof(Spectrum));
    checkmem(spec, __FILE__, __LINE__);
    spec->nSamples = nSamples;
    spec->nObs = 0;

    spec->v = malloc((nSamples-1) * sizeof(spec->v[0]));
    checkmem(spec->v);
    return spec;
}

void        Spectrum_free(Spectrum * spec) {
    free(spec->v);
    free(spec);
}

void        Spectrum_print(Spectrum * spec, FILE * ofp) {
    fprintf(ofp, "%2s %8.6lf\n", "i", "spec[i]");
    for(unsigned i=0; i < spec->nSamples; ++i)
        fprintf(ofp, "%2u %8.6lf\n", i+1, spec->v[i]/((double) spec->nObs));
}

/**
 * Sum of squared differences between two spectra
 */
double      Spectrum_sqrDiff(Spectrum *spec1, Spectrum *spec1) {
    if(spec1->nSamples != spec2->nSamples)
        eprintf("%s:%s:%d: spectra must have equal sample sizes\n",
                __FILE__, __func__, __LINE__);
    double ss = 0.0;
    for(unsigned i=1; i < spec1->nSamples; ++i) {
        double diff = Spectrum_value(spec1, i) - Spectrum_value(spec2, i);
        ss += diff*diff;
    }
    return ss;
}

