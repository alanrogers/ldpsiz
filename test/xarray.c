/**
 * @file xarray.c
 * @author Alan R. Rogers
 * @brief Test array.c
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "array.h"
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {


    if(argc != 1)
        eprintf("usage: xarray\n");

    unsigned long i, nitems = 10;

    UIntArray *ua = UIntArray_new(nitems);
    UIntArray *ua2 = UIntArray_new(nitems);
    assert(nitems == UIntArray_dim(ua));
    for(i=0; i<nitems; ++i)
        UIntArray_set(ua, i, (unsigned) i);
    for(i=0; i<nitems; ++i)
        assert(((unsigned) i) == UIntArray_get(ua, i));
    UIntArray_copy(ua2, ua);
    UIntArray *ua3 = UIntArray_dup(ua2);
    unsigned *uptr = UIntArray_ptr(ua);
    for(i=0; i<nitems; ++i) {
        assert(UIntArray_get(ua2, i) == UIntArray_get(ua, i));
        assert(UIntArray_get(ua3, i) == UIntArray_get(ua, i));
        assert(uptr[i] == UIntArray_get(ua, i));
    }
    UIntArray_free(ua);
    UIntArray_free(ua2);
    UIntArray_free(ua3);
    ua = ua2 = ua3 = NULL;
    unitTstResult("UIntArray", "OK");

    ULIntArray *ula = ULIntArray_new(nitems);
    ULIntArray *ula2 = ULIntArray_new(nitems);
    assert(nitems == ULIntArray_dim(ula));
    for(i=0; i<nitems; ++i)
        ULIntArray_set(ula, i, (unsigned long) i);
    for(i=0; i<nitems; ++i)
        assert(((unsigned long) i) == ULIntArray_get(ula, i));
    ULIntArray_copy(ula2, ula);
    ULIntArray *ula3 = ULIntArray_dup(ula2);
    unsigned long *ulptr = ULIntArray_ptr(ula);
    for(i=0; i<nitems; ++i) {
        assert(ULIntArray_get(ula2, i) == ULIntArray_get(ula, i));
        assert(ULIntArray_get(ula3, i) == ULIntArray_get(ula, i));
        assert(ulptr[i] == ULIntArray_get(ula, i));
    }
    ULIntArray_free(ula);
    ULIntArray_free(ula2);
    ULIntArray_free(ula3);
    ula = ula2 = ula3 = NULL;
    unitTstResult("ULIntArray", "OK");

    IntArray *ia = IntArray_new(nitems);
    IntArray *ia2 = IntArray_new(nitems);
    assert(nitems == IntArray_dim(ia));
    for(i=0; i<nitems; ++i)
        IntArray_set(ia, i, (int) i);
    for(i=0; i<nitems; ++i)
        assert(((int) i) == IntArray_get(ia, i));
    IntArray_copy(ia2, ia);
    IntArray *ia3 = IntArray_dup(ia2);
    int *iptr = IntArray_ptr(ia);
    for(i=0; i<nitems; ++i) {
        assert(IntArray_get(ia2, i) == IntArray_get(ia, i));
        assert(IntArray_get(ia3, i) == IntArray_get(ia, i));
        assert(iptr[i] == IntArray_get(ia, i));
    }
    IntArray_free(ia);
    IntArray_free(ia2);
    IntArray_free(ia3);
    ia = ia2 = ia3 = NULL;
    unitTstResult("IntArray", "OK");

    DblArray *da = DblArray_new(nitems);
    DblArray *da2 = DblArray_new(nitems);
    assert(nitems == DblArray_dim(da));
    for(i=0; i<nitems; ++i)
        DblArray_set(da, i, (double) i);
    for(i=0; i<nitems; ++i)
        assert(((double) i) == DblArray_get(da, i));
    DblArray_copy(da2, da);
    DblArray *da3 = DblArray_dup(da2);
    double *dptr = DblArray_ptr(da);
    for(i=0; i<nitems; ++i) {
        assert(DblArray_get(da2, i) == DblArray_get(da, i));
        assert(DblArray_get(da3, i) == DblArray_get(da, i));
        assert(dptr[i] == DblArray_get(da, i));
    }
    DblArray_free(da);
    DblArray_free(da2);
    DblArray_free(da3);
    da = da2 = da3 = NULL;
    unitTstResult("DblArray", "OK");

    return 0;
}
