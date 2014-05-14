/**
 * @file locref.h
 * @brief Header for locref.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef ARR_LOCREF
#define ARR_LOCREF
#define LOCREFMAX 100

#include <stddef.h>
#include "typedefs.h"

/**
 * LocRef facilitates using a single memory block to allocate a
 * structure and all the arrays within it. This trick is called the
 * "struct hack" in C programming.
 */
struct LocRef {
    int         n, curr;
    size_t      start[LOCREFMAX], end[LOCREFMAX];
};

void        LocRef_init(LocRef * lr);
void        LocRef_add(LocRef * lr, size_t size);
size_t      LocRef_size(LocRef * lr);
void        LocRef_setNext(LocRef * lr, int i);
void       *LocRef_next(LocRef * lr, void *base);

#endif
