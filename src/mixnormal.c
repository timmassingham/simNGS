#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "utility.h"
#include "normal.h"
#include "random.h"
#include "mixnormal.h"


void free_NormMixParam(NormMixParam param){
	if(NULL==param){ return; }
	free(param->mean);
	free(param->sd);
	free(param->prob);
	free(param);
}

NormMixParam new_NormMixParam( const unsigned int nmix){
	NormMixParam param = calloc(1,sizeof(*param));
	if(NULL==param){ return NULL;}
	param->nmix = nmix;
	param->mean = calloc(nmix,sizeof(real_t));
	param->sd = calloc(nmix,sizeof(real_t));
	param->prob = calloc(nmix,sizeof(real_t));

	if (NULL==param->mean || NULL==param->sd || NULL==param->prob ){
		free_NormMixParam(param);
		param = NULL;
	}
	return param;
}

NormMixParam copy_NormMixParam(const NormMixParam param){
	if(NULL==param){ return NULL;}
	NormMixParam newparam = new_NormMixParam(param->nmix);
	if(NULL==newparam){ return NULL; }
	for ( int i=0 ; i<param->nmix ; i++){
		newparam->prob[i] = param->prob[i];
		newparam->mean[i] = param->mean[i];
		newparam->sd[i] = param->sd[i];
	}
	return newparam;
}

NormMixParam initialise_mixture(const real_t * x, const unsigned int n, const unsigned int nmix){
	if(NULL==x){ return NULL;}
	NormMixParam param = new_NormMixParam(nmix);
	if(NULL==param){ return NULL; }

	real_t mean = 0.0;
	real_t var = 0.0;
	for ( int i=0 ; i<n ; i++){
		mean += x[i];
		var += x[i]*x[i];
	}
	mean /= n;
	var = var/n - mean*mean;
	real_t sd = sqrt(var);

	for ( int j=0 ; j<nmix ; j++){
		param->prob[j] = 1.0/nmix;
		param->mean[j] = mean + (2.0-4.0*((j+1.0)/(nmix+1.0))) * sd;
		param->sd[j] = sd/nmix;
	}

	return param;
}


real_t emstep_mixnormal(const real_t * x, const unsigned int n, NormMixParam param){
	if(NULL==param){ return NAN; }
	if(NULL==x){ return NAN; }
	const unsigned int nmix = param->nmix;
	real_t m_mean[nmix];
	real_t m_var[nmix];
	real_t m_we[nmix];
	real_t logprob[nmix];
	real_t logsd[nmix];
	for ( int j=0 ; j<nmix ; j++){
		m_mean[j] = m_var[j] = m_we[j] = 0.0;
		logprob[j] = log(param->prob[j]);
		logsd[j] = log(param->sd[j]);
	}
	real_t loglike = 0.0;
	for ( int i=0 ; i<n ; i++){
		// E-step, conditional probability
		real_t cprob[nmix];
		real_t tot = 0.0;
		real_t maxp = -HUGE_VAL;
		for ( int j=0 ; j<nmix ; j++){
			real_t tstat = (x[i]-param->mean[j])/param->sd[j];
			cprob[j] = logprob[j] -0.5*tstat*tstat - logsd[j];
			if(cprob[j]>maxp){ maxp = cprob[j];}
		}
		for ( int j=0 ; j<nmix ; j++){
			cprob[j] = exp(cprob[j]-maxp);
			tot += cprob[j];
		}
		loglike += log(tot) + maxp;
		for ( int j=0 ; j<nmix ; j++){
			cprob[j] /= tot;
		}

		// M-step, accumulate for estimates of mixture mean and variance
		for ( int j=0 ; j<nmix ; j++){
			m_mean[j] += cprob[j] * x[i];
			m_var[j] += cprob[j] * x[i] * x[i];
			m_we[j] += cprob[j];
		}
	}
	
	for ( int j=0 ; j<nmix ; j++){
		param->mean[j] = m_mean[j] / m_we[j];
		param->sd[j] = sqrt(m_var[j]/m_we[j] - param->mean[j]*param->mean[j]);
		param->prob[j] = m_we[j]/n;
	}

	return loglike - 0.5 * log(2.0*M_PI) * n;
}

void show_NormMixParam( FILE * fp, const NormMixParam param){
	fprintf(fp,"Normal mixture, %u mixtures.\n",param->nmix);
	for ( int i=0 ; i<param->nmix ; i++){
		fprintf(fp,"%3d: %6.4f %6.4f %6.4f\n",i+1,param->prob[i],param->mean[i],param->sd[i]);
	}
}


NormMixParam fit_mixnormal(const real_t * x, const unsigned int n, const unsigned int nmix, const unsigned niter){
	if(NULL==x){ return NULL;}
	NormMixParam param = initialise_mixture(x,n,nmix);
	if(NULL==param){ return NULL; }
	for ( int i=0 ; i<niter ; i++){
		emstep_mixnormal(x,n,param);
	}
	return param;
}


// Could be improved from a rounding perspective
real_t dmixnorm(const real_t x, const NormMixParam param, const bool logd){
	if(NULL==param){ return NAN;}
	real_t res = 0.0;
	for ( int i=0 ; i<param->nmix ; i++){
		res += param->prob[i] * dnorm(x,param->mean[i],param->sd[i],false);
	}
	return logd?log(res):res;
}

real_t rmixnorm(const NormMixParam param){
	if(NULL==param){ return NAN;}
	int c = rchoose(param->prob,param->nmix);
	return rnorm(param->mean[c],param->sd[c]);
}

// Could be improved from a rounding perspective
real_t pmixnorm(const real_t x, const NormMixParam param, const bool tail, const bool logp){
	if(NULL==param){ return NAN;}
	real_t res = 0.0;
	for ( int i=0 ; i<param->nmix ; i++){
		res += param->prob[i] * pnorm(x,param->mean[i],param->sd[i],tail,false);
	}
	return logp?log(res):res;
}


// Newton's method
real_t qmixnorm(const real_t p, const NormMixParam param, const bool tail, const bool logp){
	const int it_max = 20;
	const real_t tol= 3e-8;
	if(NULL==param){ return NAN;}

	// Initialise
	real_t ap = logp?exp(p):p;
	// Crude fit to logistic to find initial x
	real_t im=0.0, im2=0.0, ivar=0.0;
	for ( int i=0 ; i<param->nmix ; i++){
		im += param->prob[i] * param->mean[i];
		im2 += param->prob[i] * param->mean[i] * param->mean[i];
		ivar += param->prob[i] * param->sd[i] * param->sd[i];
	}
	real_t var = ivar + (im2-im*im);

	real_t sd = sqrt(3.0*var)/M_PI;
	real_t x = im + sd * (log(ap)-log1p(-ap));

	
	// Iteration
	for ( int i=0 ; i<it_max ; i++){
		real_t delta = (ap-pmixnorm(x,param,tail,false))/dmixnorm(x,param,false);
		x += delta;
		if(fabs(delta)/(x+3e-8) < tol){ break;}
	}
	return x;
}


#ifdef TEST

int main ( int argc, char * argv[]){
	if(argc!=4){
		errx(EXIT_FAILURE,"Usage: prog nmix niter file");
	}
	const unsigned nmix = atoi(argv[1]);
	const unsigned niter = atoi(argv[2]);
	FILE * fp = fopen(argv[3],"r");
	if(NULL==fp){ errx(EXIT_FAILURE,"Failed to open %s\n",argv[3]);}

	unsigned int nelt = 0;
	fscanf(fp,"%u",&nelt);
	real_t * x = calloc(nelt,sizeof(real_t));
	for ( int i=0 ; i<nelt ; i++){
		fscanf(fp,real_format_str,&x[i]);
	}
	fclose(fp);
        
	NormMixParam param = initialise_mixture(x,nelt,nmix);
	for ( int i=0 ; i<niter ; i++){
		show_NormMixParam(stdout,param);
		real_t l = emstep_mixnormal(x,nelt,param);
		fprintf(stdout,"loglike = %f\n",l);
	}
	show_NormMixParam(stdout,param);

	int skip = nelt/10;
	for ( int i=0 ; i<nelt ; i+=skip){
		real_t p = pmixnorm(x[i],param,false,false);
		fprintf(stdout,"x=%f\tp=%f\tx'=%f\n",x[i],p,qmixnorm(p,param,false,false));
	}

        return EXIT_SUCCESS;
}
#endif /* TEST */

