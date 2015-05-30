/**
   @file espectrum.c
   @brief Expected site frequency spectrum under piecewise-constant population size.
   @internal
 
   @copyright Copyright (c) 2014, Alan R. Rogers
   This file is released under the Internet Systems Consortium
   License, which can be found in file "LICENSE".
 
   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
**/

#include "espectrum.h"
#include "polya.h"
#include "matcoal.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>

/// Site frequency spectrum
struct ESpectrum {
	unsigned nSamples;
	double *spec;
};

/**
   Class ESpectrum: Expected site frequency spectrum under random
   mating with a given history of population size.

   @author Alan Rogers
   @date 27 Nov 2014
**/
/// Constructor
/// Polya_prob(polya,i,k) = probability that a mutation occurring in
/// the coalescent interval containing k lineages will have
/// i descendants in the modern sample.
ESpectrum *ESpectrum_new(unsigned nSamples, PopHist *ph,
                         const Polya *polya, double errTol) {
	ESpectrum *spectrum = malloc(sizeof(ESpectrum));
	checkmem(spectrum, __FILE__, __LINE__);

	spectrum->nSamples = nSamples;
	spectrum->spec = malloc((nSamples-1) * sizeof(spectrum->spec[0]));
	checkmem(spectrum->spec, __FILE__, __LINE__);

	// m[i] is expected length of coalescent interval
	// containing i+1 lineages.
	double m[nSamples];
	MatCoal_integrate(nSamples, m, ph, errTol);

	unsigned i, j, k;
	double sum=0.0;
    unsigned remainder = (nSamples-1) % 5;  // for unrolled loop
	for(i=1; i < nSamples; ++i) {
		spectrum->spec[i-1] = 0;
		// spectrum->spec[i-1] is probability that a polymorphic
		// site will have i copies of the mutant allele.
#if 0
        // Normal loop
		for(j=0; j<nSamples-1; ++j) {
            k = j+2;
			spectrum->spec[i-1] += k * m[k-1] * Polya_prob(polya, i, k);
        }
#else
        // Unrolled
        for(j=0; j<remainder; ++j) {
            k = j+2;
			spectrum->spec[i-1] += k * m[k-1] * Polya_prob(polya, i, k);
        }
        for(j=remainder; j < nSamples-1; j+=5) {
            k = j+2;
			spectrum->spec[i-1] +=
                k * m[k-1] * Polya_prob(polya, i, k)
                + (k+1) * m[k] * Polya_prob(polya, i, k+1)
                + (k+2) * m[k+1] * Polya_prob(polya, i, k+2)
                + (k+3) * m[k+2] * Polya_prob(polya, i, k+3)
                + (k+4) * m[k+3] * Polya_prob(polya, i, k+4);
        }
#endif        
		sum += spectrum->spec[i-1];
	}

	// normalize
	for(i=0; i < -1 + spectrum->nSamples; ++i)
		spectrum->spec[i] /= sum;

	return spectrum;
}

/// Destructor
void ESpectrum_free(ESpectrum *spectrum) {
	assert(spectrum);
	free(spectrum->spec);
	free(spectrum);
}

/// Return probability that mutant occurs i times in sample
/// float_type
double ESpectrum_unfolded(ESpectrum *spectrum, unsigned i) {
	assert(i > 0);
	assert(i < spectrum->nSamples);
	return spectrum->spec[i-1];
}

/// Return probability that mutant occurs either i or n-i
/// times in sample.
///
/// @param i an integer in the range 1..(_n+1)/2
double ESpectrum_folded(ESpectrum *spectrum, unsigned i) {
	assert(i > 0);
	assert(i < spectrum->nSamples);
	unsigned two_i = 2*i;
	if(two_i == spectrum->nSamples) // at crease, where spectrum is folded
		return spectrum->spec[i-1];
	if(two_i > spectrum->nSamples)
		return 0.0;
	return(spectrum->spec[i-1]
		   + spectrum->spec[spectrum->nSamples-i-1]);
}

/// Return the number of genes in the sample.
unsigned ESpectrum_nSamples(const ESpectrum *spectrum) {
	return spectrum->nSamples;
}

