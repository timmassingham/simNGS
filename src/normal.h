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

#ifndef _NORMAL_H
#define _NORMAL_H

#include "utility.h"
#include <stdint.h>
#include "matrix.h"

real_t rstdnorm( void );
real_t dstdnorm( const real_t x, const bool logd);
real_t pstdnorm( const real_t q, const bool tail, const bool logd);
real_t qstdnorm(const real_t p, const bool tail, const bool logd);

real_t rnorm( const real_t mean, const real_t sd );
real_t dnorm(  real_t x,  real_t m,  real_t sd,  bool logd);
real_t pnorm( const real_t q, const real_t m, const real_t sd, const bool tail, const bool logd);
real_t qnorm(const real_t p, const real_t m, const real_t sd, const bool tail, const bool logd);

MAT rmultinorm( const MAT mean, const MAT At, const uint32_t n, MAT z);
real_t dmultinorm( const MAT x, const MAT mean, const MAT At, const uint32_t n, const bool logd);

real_t rlognorm ( const real_t logmean, const real_t logsd);

#endif

