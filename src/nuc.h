/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the simNGS software for simulating likelihoods
 *  for next-generation sequencing machines.
 *
 *  simNGS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  simNGS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with simNGS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NUC_H
#define _NUC_H

#include <stdio.h>
#include <stdint.h>

typedef char NUC;
#define NBASE       4

#define NUC_AMBIG   4
#define NUC_A       0
#define NUC_C       1
#define NUC_G       2
#define NUC_T       3

typedef char PHREDCHAR;

void show_NUC(FILE * fp, const NUC nuc);
    
NUC nuc_from_char( const char c);
char char_from_nuc(const NUC nuc);
NUC complement(const NUC nuc);
NUC * reverse_complement(const NUC * nuc, const uint32_t len);
PHREDCHAR phredchar_from_char( const char c);

#endif

