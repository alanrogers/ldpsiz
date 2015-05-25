#include "array.h"
#include <string.h>

///////////////////////////////////////////////////////

UIntArray *UIntArray_new(unsigned long nItems) {
    if(nItems == 0)
        eprintf("%s:%s:%d: can't allocate an empty array",
                __FILE__, __func__, __LINE__);
    size_t size = (nItems-1)*sizeof(unsigned) + sizeof(UIntArray);
    UIntArray *new = malloc(size);
    checkmem(new, __FILE__, __LINE__);
    new->nItems = nItems;
    return new;
}

void UIntArray_free(UIntArray *a) {
    free(a);
}

void        UIntArray_copy(UIntArray * to, const UIntArray *from) {
#ifndef NDEBUG
    if(to->nItems != from->nItems)
        eprintf("%s:%s:%d: arrays are of unequal size",
                __FILE__, __func__, __LINE__);
#endif
    memcpy(to->v, from->v, to->nItems * sizeof(to->v[0]));
}

///////////////////////////////////////////////////////

ULIntArray *ULIntArray_new(unsigned long nItems) {
    if(nItems == 0)
        eprintf("%s:%s:%d: can't allocate an empty array",
                __FILE__, __func__, __LINE__);
    size_t size = (nItems-1)*sizeof(unsigned long) + sizeof(ULIntArray);
    ULIntArray *new = malloc(size);
    checkmem(new, __FILE__, __LINE__);
    new->nItems = nItems;
    return new;
}

void ULIntArray_free(ULIntArray *a) {
    free(a);
}

void        ULIntArray_copy(ULIntArray * to, const ULIntArray *from) {
#ifndef NDEBUG
    if(to->nItems != from->nItems)
        eprintf("%s:%s:%d: arrays are of unequal size",
                __FILE__, __func__, __LINE__);
#endif
    memcpy(to->v, from->v, to->nItems * sizeof(to->v[0]));
}

///////////////////////////////////////////////////////

IntArray *IntArray_new(unsigned long nItems) {
    if(nItems == 0)
        eprintf("%s:%s:%d: can't allocate an empty array",
                __FILE__, __func__, __LINE__);
    size_t size = (nItems-1)*sizeof(int) + sizeof(IntArray);
    IntArray *new = malloc(size);
    checkmem(new, __FILE__, __LINE__);
    new->nItems = nItems;
    return new;
}

void IntArray_free(IntArray *a) {
    free(a);
}

void        IntArray_copy(IntArray * to, const IntArray *from) {
#ifndef NDEBUG
    if(to->nItems != from->nItems)
        eprintf("%s:%s:%d: arrays are of unequal size",
                __FILE__, __func__, __LINE__);
#endif
    memcpy(to->v, from->v, to->nItems * sizeof(to->v[0]));
}

///////////////////////////////////////////////////////

DblArray *DblArray_new(unsigned long nItems) {
    if(nItems == 0)
        eprintf("%s:%s:%d: can't allocate an empty array",
                __FILE__, __func__, __LINE__);
    size_t size = (nItems-1)*sizeof(double) + sizeof(DblArray);
    DblArray *new = malloc(size);
    checkmem(new, __FILE__, __LINE__);
    new->nItems = nItems;
    return new;
}

void DblArray_free(DblArray *a) {
    free(a);
}

void        DblArray_copy(DblArray * to, const DblArray *from) {
#ifndef NDEBUG
    if(to->nItems != from->nItems)
        eprintf("%s:%s:%d: arrays are of unequal size",
                __FILE__, __func__, __LINE__);
#endif
    memcpy(to->v, from->v, to->nItems * sizeof(to->v[0]));
}
