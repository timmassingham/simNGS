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
 
#ifndef _KUMARASWAMY_H
#define _KUMARASWAMY_H

#include <stdbool.h>
#include "utility.h"

real_t pkumaraswamy(  real_t x,  real_t shape1,  real_t shape2,  bool tail,  bool logp);
real_t qkumaraswamy(  real_t p,  real_t shape1,  real_t shape2,  bool tail,  bool logp);
real_t dkumaraswamy(  real_t x,  real_t shape1,  real_t shape2,  bool logd );
real_t rkumaraswamy(  real_t shape1,  real_t shape2 );


#endif /* _KUMARASWAMY_H */

