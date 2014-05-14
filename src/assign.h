/**
 * @file assign.h
 * @author Alan R. Rogers
 * @brief Header for functions defined in assign.c
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_ASSIGNMENT_H
#define LDPSIZ_ASSIGNMENT_H

#include "typedefs.h"
#include <stdio.h>
#define MANDATORY 1

Assignment *Assignment_new(Assignment * a, const char *key,
                           const char *value);
void        Assignment_free(Assignment * a);
const char *Assignment_value(const Assignment * a, const char *key);
void        Assignment_print(const Assignment * a, FILE * fp);
int         Assignment_setDbl(const Assignment * a, const char *key,
                              double *ptr, int mandatory);
int         Assignment_setLong(const Assignment * a, const char *key,
                               long *ptr, int mandatory);
int         Assignment_setInt(const Assignment * a, const char *key, int *ptr,
                              int mandatory);
unsigned    Assignment_setUnsignedInt(const Assignment * a, const char *key,
                                      unsigned *ptr, int mandatory);
int         Assignment_setString(const Assignment * a, const char *key,
                                 char *ptr, int size, int mandatory);
Assignment *Assignment_readObsld(FILE * ifp);

#endif
