/**
 * @file ini.h
 * @author Alan R. Rogers
 * @brief Header for ini.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LDPSIZ_INI_H
#define LDPSIZ_INI_H

#include "typedefs.h"
#include <stdio.h>
#define MANDATORY 1

Ini        *Ini_new(const char *ifname);
void        Ini_print(Ini * ini, FILE * ofp);
void        Ini_free(Ini * ini);

double      Ini_twoN0(const Ini * ini);
int         Ini_setDbl(const Ini * ini, const char *key, double *ptr,
                       int mandatory);
int         Ini_setLong(const Ini * ini, const char *key, long *ptr,
                        int mandatory);
int         Ini_setInt(const Ini * ini, const char *key, int *ptr,
                       int mandatory);
unsigned    Ini_setUnsignedInt(const Ini * ini, const char *key,
                               unsigned *ptr, int mandatory);
int         Ini_setString(const Ini * ini, const char *key, char *ptr,
                          int size, int mandatory);
int         Ini_setEpochLink(const Ini * ini, EpochLink ** ptr,
                             int mandatory);

#endif
