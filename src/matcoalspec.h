#ifndef MATCOALSPEC_INCLUDED
#  define MATCOALSPEC_INCLUDED

#  include "typedefs.h"
#  include <mpfr.h>

struct MpfrVec {
    unsigned    dim;
    mpfr_t     *x;
};

MpfrVec    *MpfrVec_new(unsigned dim, long double x[dim]);
void        MpfrVec_free(MpfrVec * self);
void        MpfrVec_get(MpfrVec *self, unsigned dim, long double x[dim]);

MatCoalSpec *MatCoalSpec_new(unsigned nSamples);
void        MatCoalSpec_print(MatCoalSpec * self);
void        MatCoalSpec_project(MatCoalSpec * self, MpfrVec * x,
                                long double v);

void UTmatXvec(unsigned dim, mpfr_t *y, mpfr_t *A, unsigned *offset, mpfr_t *x);

#endif
