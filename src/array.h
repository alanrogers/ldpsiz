#ifndef LDPSIZ_ARRAY_H
#  define LDPSIZ_ARRAY_H

#  include "misc.h"
#  include "typedefs.h"

struct UIntArray {
    unsigned long nItems;
    unsigned    v[1];
};

struct IntArray {
    unsigned long nItems;
    int         v[1];
};

struct ULIntArray {
    unsigned long nItems;
    unsigned long v[1];
};

struct DblArray {
    unsigned long nItems;
    double        v[1];
};

UIntArray  *UIntArray_new(unsigned long nItems);
void        UIntArray_free(UIntArray * a);
void        UIntArray_copy(UIntArray * to, const UIntArray *from);
static inline unsigned long UIntArray_dim(const UIntArray *a);
static inline unsigned *UIntArray_ptr(UIntArray *a);
static inline unsigned UIntArray_get(const UIntArray * a, unsigned long i);
static inline void UIntArray_set(UIntArray * a, unsigned long i,
                                 unsigned value);

ULIntArray  *ULIntArray_new(unsigned long nItems);
void        ULIntArray_free(ULIntArray * a);
void        ULIntArray_copy(ULIntArray * to, const ULIntArray *from);
static inline unsigned long ULIntArray_dim(const ULIntArray *a);
static inline unsigned long *ULIntArray_ptr(ULIntArray *a);
static inline unsigned long ULIntArray_get(const ULIntArray * a, unsigned long i);
static inline void ULIntArray_set(ULIntArray * a, unsigned long i,
                                 unsigned long value);

IntArray   *IntArray_new(unsigned long nItems);
void        IntArray_free(IntArray * a);
void        IntArray_copy(IntArray * to, const IntArray *from);
static inline unsigned long IntArray_dim(const IntArray *a);
static inline int *IntArray_ptr(IntArray *a);
static inline int IntArray_get(const IntArray * a, unsigned long i);
static inline void IntArray_set(IntArray * a, unsigned long i,
                                int value);

DblArray   *DblArray_new(unsigned long nItems);
void        DblArray_free(DblArray * a);
void        DblArray_copy(DblArray * to, const DblArray *from);
static inline unsigned long DblArray_dim(const DblArray *a);
static inline double *DblArray_ptr(DblArray *a);
static inline double DblArray_get(const DblArray * a, unsigned long i);
static inline void DblArray_set(DblArray * a, unsigned long i,
                                double value);

///////////////////////////////////////////////////////

static inline unsigned long UIntArray_dim(const UIntArray *a) {
    return a->nItems;
}

static inline unsigned *UIntArray_ptr(UIntArray *a) {
    return a->v;
}

static inline unsigned UIntArray_get(const UIntArray * a, unsigned long i) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: reading out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    return a->v[i];
}

static inline void UIntArray_set(UIntArray * a, unsigned long i,
                                 unsigned value) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: writing out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    a->v[i] = value;
}

/////////////////////////////////////////////////

static inline unsigned long ULIntArray_dim(const ULIntArray *a) {
    return a->nItems;
}

static inline unsigned long *ULIntArray_ptr(ULIntArray *a) {
    return a->v;
}

static inline unsigned long ULIntArray_get(const ULIntArray * a, unsigned long i) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: reading out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    return a->v[i];
}

static inline void ULIntArray_set(ULIntArray * a, unsigned long i,
                                 unsigned long value) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: writing out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    a->v[i] = value;
}

///////////////////////////////////////////////////////

static inline unsigned long IntArray_dim(const IntArray *a) {
    return a->nItems;
}

static inline int *IntArray_ptr(IntArray *a) {
    return a->v;
}

static inline int IntArray_get(const IntArray * a, unsigned long i) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: reading out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    return a->v[i];
}
static inline void IntArray_set(IntArray * a, unsigned long i,
                                int value) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: writing out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    a->v[i] = value;
}

///////////////////////////////////////////////////////

static inline unsigned long DblArray_dim(const DblArray *a) {
    return a->nItems;
}

static inline double *DblArray_ptr(DblArray *a) {
    return a->v;
}

static inline double DblArray_get(const DblArray * a, unsigned long i) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: reading out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    return a->v[i];
}

static inline void DblArray_set(DblArray * a, unsigned long i,
                                 double value) {
#  ifndef NDEBUG
    if(i >= a->nItems)
        eprintf("%s:%s:%d: writing out of bounds",
                __FILE__, __func__, __LINE__);
#  endif
    a->v[i] = value;
}

#endif
