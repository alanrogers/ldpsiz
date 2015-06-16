/**
   @file tfespectrum.c
   @brief Truncated, folded expected site frequency spectrum under
   piecewise-constant population size.  
   @internal
 
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "tfespectrum.h"
#include "espectrum.h"
#include "spectab.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/// Truncated folded site frequency spectrum
struct TFESpectrum {
    unsigned truncSFS; // number of entries to truncate
    unsigned dim;
    double *tfespec;
    ESpectrum *espec;
};

TFESpectrum  *TFESpectrum_new(unsigned nSamples, unsigned truncSFS, PopHist *ph,
                            const Polya *polya, double errTol) {
    TFESpectrum *new = malloc(sizeof(*new));
    checkmem(new, __FILE__, __LINE__);

    new->truncSFS = truncSFS;
    new->dim = specdim(nSamples, true);

    new->tfespec = malloc(nSamples * sizeof(new->tfespec[0]));
    checkmem(new->tfespec, __FILE__, __LINE__);

    new->espec = ESpectrum_new(nSamples, ph, polya, errTol);
    checkmem(new->espec, __FILE__, __LINE__);

    unsigned i;
    double sum=0.0;
    assert(truncSFS < new->dim);
    for(i = truncSFS; i < new->dim; ++i) {
        new->tfespec[i] = ESpectrum_folded(new->espec, i+1);
        sum += new->tfespec[i];
    }
    for(i = truncSFS; i < new->dim; ++i)
        new->tfespec[i] /= sum;
    return new;
}

void        TFESpectrum_free(TFESpectrum * self) {
    assert(self);
    ESpectrum_free(self->espec);
    free(self->tfespec);
    free(self);
}

double      TFESpectrum_atNdx(TFESpectrum * self, unsigned i) {
	assert(i > 0);
    assert(i < 1 + self->dim);
    return self->tfespec[i-1];
}

unsigned    TFESpectrum_nSamples(const TFESpectrum * self) {
    return ESpectrum_nSamples(self->espec);
}

/// Return sum of squared differences between self and an unnormalized
/// spectrum.
double TFESpectrum_diff(const TFESpectrum *self, unsigned dim, double s[dim]) {
    assert(self);
    assert(dim == self->dim);
    unsigned i;
    double ssq=0.0, sum=0.0;
    for(i = self->truncSFS; i < dim; ++i)
        sum += s[i];
    for(i = self->truncSFS; i < self->dim; ++i) {
        double diff = self->tfespec[i] - s[i]/sum;
        ssq += diff*diff;
    }
    return ssq;
}

/// Return Kullback-Leibler divergence between self and an estimated spectrum.
double TFESpectrum_KLdiverg(const TFESpectrum *self, unsigned dim, double s[dim]) {
    assert(self);
    assert(dim == self->dim);
    unsigned i;
    double kl=0.0, sum=0.0;
    for(i = self->truncSFS; i < dim; ++i)
        sum += s[i];
    for(i = self->truncSFS; i < self->dim; ++i) {
        double q = s[i]/sum;
        double p =  self->tfespec[i];
        kl += p*log(p/q);
    }
    return kl;
}

