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

#include <stdbool.h>
#include "utility.h"

real_t pweibull(  real_t x,  real_t shape,  real_t scale,  bool tail,  bool logp);
real_t qweibull(  real_t p,  real_t shape,  real_t scale,  bool tail,  bool logp);
real_t dweibull(  real_t x,  real_t shape,  real_t scale,  bool logd );
real_t rweibull(  real_t shape,  real_t scale );

#endif
