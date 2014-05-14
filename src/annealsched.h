/* 
 * Copyright (C) 2014 Alan Rogers
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
#ifndef __ARR_ANNEALSCHED__
#define __ARR_ANNEALSCHED__

typedef struct AnnealSched AnnealSched;

#include <stdio.h>

AnnealSched *AnnealSched_alloc(int nT, double initT, double decay);
int          AnnealSched_size(const AnnealSched *sched);
AnnealSched *AnnealSched_copy(const AnnealSched *old);
int          AnnealSched_cmp(const AnnealSched *s1, const AnnealSched *s2);
void         AnnealSched_reset(AnnealSched *s);
double       AnnealSched_next(AnnealSched * s);
void         AnnealSched_free(AnnealSched * s);
void         AnnealSched_print(AnnealSched *s, FILE *fp);
#endif /* __ARR_ANNEALSCHED__ */
