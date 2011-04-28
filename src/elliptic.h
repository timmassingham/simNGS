/*
 *  Copyright (C) 2011 by Tim Massingham, European Bioinformatics Institute
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

#ifndef _ELLIPTIC_H
#define _ELLIPTIC_H

#include "matrix.h"
#include "utility.h"

// Radius functions
real_t lognormal_radius( int n );
real_t normal_radius( int n);

MAT relliptic ( const MAT mean, const MAT L, real_t (*randomradius)(int), const uint32_t n, MAT z);

#endif

