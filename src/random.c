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

#include <stdint.h>
#include <assert.h>
#include "utility.h"
#include "random.h"
#include <math.h>

uint32_t rchoose( const real_t * p, const uint32_t n){
    real_t x = runif();
    uint32_t i=0;
    while( i<n && x>p[i]){ x-=p[i]; i++; }
    return (i<n)?i:(n-1); // Incase of numeric problems
}

/** Random gamma variable when shape>1
  * See Knuth V2, 3.41.
  */
real_t rgamma1( const real_t shape){
        assert(shape>1);
        real_t v,x,y;
        const real_t shm1 = shape - 1.0;
        const real_t sh2m1 = sqrt(2.0*shape-1.0);

start:
        y = tan(M_PI*genrand_real3());
        x = sh2m1*y + shm1;

        if(x<=0.0){ goto start;}
        
        v = genrand_real3();
        if(v>(1.0+y*y)*exp(shm1*log(x/shm1)-sh2m1*y)){ goto start;}
        return x;
}

/** Random gamma variable when shape<1
  * See Kundu and Gupta 2007
  * "A convenient way of generating gamma random variables using generalized exponential distribution"
  */
real_t rgamma2 ( const real_t shape){
        assert(shape>0 && shape<1);
        real_t u,v,x;
        const real_t d = 1.0334 - 0.0766*exp(2.2942*shape);
        const real_t a = exp2(shape) * pow(-expm1(-d/2),shape);
        const real_t pdsh = pow(d,shape-1.0);
        const real_t b = shape * pdsh * exp(-d);
        const real_t c = a+b;
        
start:
        u = genrand_real3();
        x = (u<=a/c) ? -2.0 * log1p(-pow(c*u,1.0/shape)/2.0) : -log( c*(1.0-u)/(shape*pdsh) );
        v = genrand_real3();
        if(x<=d){
                const real_t p = pow(x,shape-1.0)*exp(-x/2.0)/( exp2(shape-1.0)*pow(-expm1(-x/2.0),shape-1.0) );
                if(v>p){ goto start; }
        } else {
                if(v>pow(d/x,1.0-shape)){ goto start;}
        }
        
        return x;
}

/* Exponential distribution with rate */
real_t rexp(const real_t r){
        return -log(genrand_real3())/r;
}


real_t rgamma(const real_t shape, const real_t scale){
        assert(shape>0.0 && scale>0.0);
        if(shape==1.0){ return scale*rexp(1.0); }
        if(shape>1.0){ return scale*rgamma1(shape); }
        return scale*rgamma2(shape);
}
real_t rchisq(const real_t df){
        assert(df>0.0);
        return rgamma(df/2.0,2.0);
}


