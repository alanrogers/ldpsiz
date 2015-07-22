#ifndef MATCOALSPEC_INCLUDED
#define MATCOALSPEC_INCLUDED

#include "typedefs.h"

MatCoalSpec *MatCoalSpec_new(unsigned nSamples, unsigned precision);
void MatCoalSpec_print(MatCoalSpec *self);

#endif
