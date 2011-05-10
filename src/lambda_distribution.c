#include <stdlib.h>
#include <err.h>
#include <math.h>
#include "lambda_distribution.h"
#include "weibull.h"

int nparameter_distribution(const char dist){
	switch(dist){
		case 0: return 0;
		case 'W': return 2;
		case 'L': return 2;
		default: errx(EXIT_FAILURE,"Unrecognised distribution \"%c\" in %s",dist,__func__);
	}
	// Never reach here
	return -1;
}

bool validate_parameters(const real_t * param , const char dist){
	int np = nparameter_distribution(dist);
	if(np>0 && NULL==param){ return false;}
	switch(dist){
		case 0:   break;
		case 'W': if(param[0]<0.){ return false; }  //shape
		          if(param[1]<0.){ return false; }  //scale
			  break;
		case 'L': if(param[1]<0.){ return false; } // scale
			  break;
		default: errx(EXIT_FAILURE,"Unrecognised distribution in %s",__func__);
	}
	return true;
}

real_t qdistribution(const real_t px, const char dist, const real_t * param, const bool tail, const bool logp){
	if(NULL==param){ return NAN;}
	switch(dist){
		case 0:   return 0;
		case 'W': return qweibull(px,param[0],param[1],tail,logp);
		case 'L': return qlogistic(px,param[0],param[1],tail,logp);
		default: errx(EXIT_FAILURE,"Unrecognised distribution in %s",__func__);
	}
	// Never reach here
	return NAN;
}

real_t qlogistic(const real_t p, const real_t loc, const real_t sc, const bool tail, const bool logp){
	real_t pq = (logp)? expm1(p+M_LN2) : (2.0*p-1.0);
	if(tail){
		pq = -pq;
	}
	return 2.0*sc*atanh(pq) + loc;
}
