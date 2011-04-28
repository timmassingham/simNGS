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
 *  GNU General Public License for more details. *
 *  You should have received a copy of the GNU General Public License
 *  along with simNGS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <tgmath.h>
#include <stdbool.h>
#include "random.h"
#include "lapack.h"
#include "normal.h"


// Normal generator, squared radius is chi-sq
real_t normal_radius( int n){
	return sqrt(rchisq((real_t)n));
}

// Squared radius is log-normal
static real_t lognormal_mean = 5.07;
static real_t lognormal_sd = 1.088;
real_t lognormal_radius( int n ){
	return sqrt( rlognorm(lognormal_mean,lognormal_sd) );
}


/* Generate multivariate normal obserations
 * L is a upper/lower Cholesky factorisation of the variance matrix as produced
 * by routine "cholesky"
 * x = mean + L * z
 * where z is independent standard normal.
 * If mean==NULL, assume mean is zero
 * If L==NULL, assume independence (variance is identity matrix).
 */
MAT relliptic ( const MAT mean, const MAT L, real_t (*randomradius)(int), const uint32_t n, MAT z){
    if ( NULL==z){
        z = new_MAT(n,1);
        validate(NULL!=z,NULL);
    }
    
    
    // Generate initial random variable uniform on n-sphere
    real_t norm = 0.;
    for ( int i=0 ; i<n ; i++){
        z->x[i] = rstdnorm();
        norm += z->x[i]*z->x[i];
    }
    norm = sqrt(norm);
    for ( int i=0 ; i<n ; i++){
        z->x[i] /= norm;
    }
    
    
    if(NULL!=L){
        #ifndef USE_BLAS
            // Form L*z via z^t L^t
            // Start at rightmost column and work back so can modify z inplace
            for ( int i=(n-1) ; i>=0 ; i--){
                z->x[i] *= L->x[i*n+i];
                for ( int j=0 ; j<i ; j++){
                    z->x[i] += z->x[j] * L->x[i*n+j];
                }
            }
       #else
           int N = n;
           trmv(LAPACK_LOWER,LAPACK_NOTRANS,LAPACK_NONUNITTRI,&N,L->x,&N,z->x,LAPACK_UNIT);
       #endif
    }

    // Multiply by radius
    real_t r = randomradius(n);
    //fprintf(stderr,"Radius = %f\n",r);
    for ( int i=0 ; i<n ; i++){
        z->x[i] *= r;
    }

    if(NULL!=mean){
        for ( int i=0 ; i<n ; i++){
            z->x[i] += mean->x[i];
        }
    }
    return z;
}

#ifdef TEST
#include <stdlib.h>
#include <stdio.h>

const real_t mean_arry[] = { -1.0, 0.0, 1.0};
real_t var_arry[] =
  { 1.0, 0.2, -0.05,
    0.2, 2.0, -0.5,
    -0.05, -0.5, 3.0
  };

const real_t chol_arry[] =
  { 1.0, 0.2, -0.05,
    0.2, 1.4, -0.35,
    -0.05, -0.35, 1.695582
  };
#define NP 3


int main(int argc, char * argv[] ){
       if(argc!=3){
               fputs("Usage: test n seed\n",stderr);
               return EXIT_FAILURE;
       }

       unsigned int n=0;
       long unsigned int seed = 0;
       sscanf(argv[1],"%u",&n);
       sscanf(argv[2],"%lu",&seed);

       init_gen_rand(seed);

       MAT mean = new_MAT_from_array(3,1,mean_arry);
       MAT chol = new_MAT_from_array(3,3,chol_arry);
       MAT f = NULL;
       for ( int i=0 ; i<n ; i++){
            f = relliptic(mean,chol,lognormal_radius,NP,f);
            for ( int j=0 ; j<NP ; j++){
                fprintf(stdout,"\t%5.4f",f->x[j]);
            }
            fputc('\n',stdout);
        }

	return EXIT_SUCCESS;
}
#endif /* TEST */


