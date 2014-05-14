/**
 * @file locref.c
 * @anchor locref
 * @brief LocRef is a class that facilitates using a single memory
 * block to allocate a structure and all the arrays within it.
 * 
 * LocRef is facilitates using a single memory block to allocate a
 * structure and all the arrays within it. This trick is called the
 * "struct hack" in C programming.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <string.h>
#include "misc.h"
#include "locref.h"

/** Fill LocRef structure with zeroes. */
void LocRef_init(LocRef * lr) {
    memset(lr, 0, sizeof(LocRef));
}

/** Add an entry of specified size bytes to the LocRef table. */
void LocRef_add(LocRef * lr, size_t bytes) {
    if(lr->n == LOCREFMAX)
        eprintf("ERR@%s:%d: LocRef structure is full. LOCREFMAX=%d",
                __FILE__, __LINE__, LOCREFMAX);

    if(lr->n == 0)
        lr->start[0] = 0;
    else
        lr->start[lr->n] = lr->end[lr->n - 1];
    lr->end[lr->n] = lr->start[lr->n] + bytes;
    lr->n += 1;
}

/** Total size of memory block to be allocated. */
size_t LocRef_size(LocRef * lr) {
    if(lr->n == 0)
        return 0;
    return lr->end[lr->n - 1];
}

/** Set the index to be returned by the next call to LocRef_next. */
void LocRef_setNext(LocRef * lr, int i) {
    lr->curr = i;
}

/** Return a pointer to the next object in LocRef */
void  *LocRef_next(LocRef * lr, void *base) {

    int         i = lr->curr;
    void       *p = (void *) (((size_t) base) + lr->start[i]);

    ++lr->curr;

    return p;
}
