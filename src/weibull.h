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

#ifndef _WEIBULL_H
#define _WEIBULL_H

#ifndef _STDBOOL_H
#include <stdbool.h>
#endif

#ifndef _UTILITY_H
#include "utility.h"
#endif

real_t pweibull( const real_t x, const real_t shape, const real_t scale, const bool tail, const bool logp);
real_t qweibull( const real_t p, const real_t shape, const real_t scale, const bool tail, const bool logp);
real_t dweibull( const real_t x, const real_t shape, const real_t scale, const bool logd );
real_t rweibull( const real_t shape, const real_t scale );

#endif
