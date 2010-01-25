#include <tgmath.h>
#include <stdbool.h>
#include "utility.h"
#include "random.h"
#include "matrix.h"


// Standard normal using Box-Muller method
real_t rstdnorm( void ){
	static real_t z;
	static bool b = false;

	if ( true==b){ b = false; return z; }

	real_t u = runif();
	real_t v = runif();

	real_t r = sqrt( -2.0 * log(u) );
	z = r * cos(2.0*M_PI*v);
	b = true;

	return r*sin(2.0*M_PI*v);
}

real_t dstdnorm( const real_t x, const bool logd){
    real_t d = -0.5*(log(M_PI*2.)+x*x);
    return (false==logd)?exp(d):d;
}
    

real_t rnormal( const real_t mean, const real_t sd ){
	if(sd<0){ return NAN;}
	return mean + sd * rstdnorm();
}

real_t dnorm( const real_t x, const real_t m, const real_t sd, const bool logd){
    real_t r = dstdnorm((x-m)/sd,logd);
    return (false==logd)? r/sd : r - log(sd);
}



/* Generate multivariate normal obserations
 * L is a upper/lower Cholesky factorisation of the variance matrix as produced
 * by routine "cholesky"
 * x = mean + L * z
 * where z is independent standard normal.
 * If mean==NULL, assume mean is zero
 * If L==NULL, assume independence (variance is identity matrix).
 */
MAT rmultinorm ( const MAT mean, const MAT L, const uint32_t n, MAT z){
    if ( NULL==z){
        z = new_MAT(n,1);
        validate(NULL!=z,NULL);
    }
    
    
    // Get initial random variables
    for ( uint32_t i=0 ; i<n ; i++){
        z->x[i] = rstdnorm();
    }
    
    if(NULL!=L){       
        // Form L*z via z^t L^t
        // Start at rightmost column and work back so can modify z inplace
        for ( int i=(n-1) ; i>=0 ; i--){
            z->x[i] *= L->x[i*n+i];
            for ( uint32_t j=0 ; j<i ; j++){
                z->x[i] += z->x[j] * L->x[i*n+j];
            }
        }
    }
    
    if(NULL!=mean){
        for ( uint32_t i=0 ; i<n ; i++){
            z->x[i] += mean->x[i];
        }
    }
    return z;
}

/* Density of multivariate normal distribution
 * Computed using the inverse of the Cholesky factorisation of the variance 
 * matrix (lower diagonal).
 */
real_t dmultinorm( const MAT x, const MAT mean, const MAT invL, const uint32_t n, const bool logd){
    validate(NULL!=x,NAN);
    
    // (x-m)^t invL * invL^t (x-m) 
    real_t sum = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        real_t rsum = 0.;
        for ( uint32_t j=0 ; j<=i ; j++){
            // Checks here are not best from performance perspective
            real_t diff = (NULL!=mean)? (x->x[j]-mean->x[j]) : x->x[j];
            rsum += (NULL!=invL) ? (diff * invL->x[i*n+j]) : diff;
        }
        sum += rsum*rsum;
    }
    
    // Determinate
    real_t det = 1.;
    if(NULL!=invL){
        for ( uint32_t i=0 ; i<n ; i++){
            det *= invL->x[i*n+i] * invL->x[i*n+i];
        }
    }
    
    return (logd==false) ? exp(-0.5*sum) * pow(2.0*M_PI,-(real_t)n/2.0) * sqrt(fabs(det)) 
                         : -0.5*sum - (n/2.0)*log(2.0*M_PI) + 0.5*log(fabs(det));
}

#ifdef TEST
#include <stdlib.h>
#include <stdio.h>

#define NP 3
const real_t mean_arry[] = { -1.0, 0.0, 1.0};

real_t var_arry[] =
  { 1.0, 0.2, -0.05,
    0.2, 2.0, -0.5,
    -0.05, -0.5, 3.0
  };

  /* (Symmetrised) Cholsky factorisation of matrix
 *  1.00  0.2 -0.05
 *  0.20  2.0 -0.50
 * -0.05 -0.5  3.00
 */

const real_t chol_arry[] =
  { 1.0, 0.2, -0.05,
    0.2, 1.4, -0.35,
    -0.05, -0.35, 1.695582
  };


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
/*	fputs("# Standard normals\n",stdout);
	for ( unsigned int i=0 ; i<n ; i++){
		fprintf(stdout,"%f\n",rstdnorm());
	}
*/
	MAT mean = new_MAT_from_array(3,1,mean_arry);
	MAT var = new_MAT_from_array(3,3,var_arry);
	MAT chol = new_MAT_from_array(3,3,chol_arry);
	
	cholesky(var);
	//fputs("Cholesky decomposition of var.\n",stderr);
	//show_MAT(stderr,var,3,3);

	MAT f = NULL;
	fputs("#Multivariate normals\n",stdout);
	for ( int i=0 ; i<n ; i++){
	    f = rmultinorm(mean,chol,NP,f);
	    //fprintf(stdout,"%3d:",i+1);
	    for ( int j=0 ; j<NP ; j++){
	        fprintf(stdout,"\t%5.4f",f->x[j]);
	    }
	    fputc('\n',stdout);
	}
/*
	fputs("Inverting Cholesky decomposition.\n",stderr);
	invert_cholesky(var);
	show_MAT(stderr,var,3,3);
	fprintf(stderr,"Density %e, log density %e\n",dmultinorm( f,mean,var,NP,false), dmultinorm( f,mean,var,NP,true));
*/
	
}
#endif

