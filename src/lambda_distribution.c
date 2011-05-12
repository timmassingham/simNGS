#include <stdlib.h>
#include <err.h>
#include <math.h>
#include <string.h>
#include "lambda_distribution.h"
#include "weibull.h"
#include "normal.h"
#include "mixnormal.h"

char * nameOfDistribution(const Distribution dist){
	if(NULL==dist){ return NULL; }
	switch(dist->key){
		case 'W': return "Weibull";
		case 'L': return "Logistic";
		case 'N': return "Normal";
		case 'M': return "Mixture of normals";
		default: errx(EXIT_FAILURE,"Unrecognised distribution in %s",__func__);
	}
	// Never reach here;
	return NULL;
}

void free_Distribution(Distribution dist){
	if(NULL==dist){ return; }
	// Special cases
	switch(dist->key){
		case 'M': free_NormMixParam((NormMixParam)dist->info);
	}
	free(dist->param);
	free(dist);
}

Distribution copy_Distribution(const Distribution dist){
	if(NULL==dist){ return NULL; }
	Distribution newdist = calloc(1,sizeof(*newdist));
	newdist->key = dist->key;
	newdist->np = dist->np;
	switch(dist->key){
		case 'M': newdist->info = (void *)copy_NormMixParam((NormMixParam)dist->info);
			  if(NULL==newdist->info){ goto cleanup; }
			  break;
		default:
			  newdist->param = calloc(dist->np,sizeof(real_t));
			  if(NULL==newdist->param){ goto cleanup; }
			  memcpy(newdist->param,dist->param,dist->np*sizeof(real_t));
	}

	return newdist;

cleanup: free_Distribution(newdist);
	 return NULL;
}

Distribution new_Distribution(const char type, const real_t * param){
	if(NULL==param){ return NULL;}
	Distribution dist = calloc(1,sizeof(*dist));
	dist->key = type;
	NormMixParam normmix = NULL;
	switch(type){
		case 'M': normmix = new_NormMixParam( (int)param[0] );
			  if(NULL==normmix){ goto cleanup;}
			  dist->info = (void *)normmix;
			  for ( int i=0,j=1 ; i<normmix->nmix ; i++){
				  normmix->prob[i] = param[j++];
				  normmix->mean[i] = param[j++];
				  normmix->sd[i] = param[j++];
			  }
			  break;
		case 'L': 
		case 'W':
		case 'N': dist->np = 2;
			  dist->param = calloc(dist->np,sizeof(real_t));
			  if(NULL==dist->param){ goto cleanup; }
			  for ( int i=0 ; i<dist->np ; i++){
				  dist->param[i] = param[i];
			  }
			  break;
		default: errx(EXIT_FAILURE,"Unrecognised distribution in %s",__func__);
	}
	return dist;

cleanup:
	free_Distribution(dist);
	return NULL;
}

void show_Distribution(FILE * fp, const Distribution dist){
	if(NULL==fp || NULL==dist){ return; }

	switch(dist->key){
		case 'M': show_NormMixParam(fp,(NormMixParam)dist->info);
			  break;
		default: fputs( nameOfDistribution(dist),fp);
			 for ( int i=0 ; i<dist->np ; i++){
				 fprintf(fp,"\t%f",dist->param[i]);
			 }
			 fputc('\n',fp);
	}
}

Distribution new_Distribution_from_fp(FILE * fp){
	if(NULL==fp){ return NULL; }
	Distribution dist = calloc(1,sizeof(*dist));
	fscanf(fp,"%c",&dist->key);
	switch(dist->key){
		case 'M': fscanf(fp,"%d",&dist->np);
			  if(dist->np<1){ goto cleanup; }
			  NormMixParam normmix = new_NormMixParam(dist->np);
			  dist->info = (void *)normmix;
			  for(int i=0 ; i<dist->np ; i++){
				  int ret = fscanf(fp,real_format_str,&normmix->prob[i]);
				  ret += fscanf(fp,real_format_str,&normmix->mean[i]);
				  ret += fscanf(fp,real_format_str,&normmix->sd[i]);
				  if(3!=ret){ goto cleanup;}
			  }
			  dist->np = 3*dist->np-1;
			  break;
		case 'W': dist->np = 2;
			  break;
		case 'L': dist->np = 2;
			  break;
		case 'N': dist->np = 2;
			  break;
	}
	if(NULL==dist->info){
		dist->param = calloc(dist->np,sizeof(real_t));
		if(NULL==dist->param){ goto cleanup; }
		for ( int i=0 ; i<dist->np ; i++){
			int ret = fscanf(fp, real_format_str , &dist->param[i]);
			if(1!=ret){ goto cleanup; }
		}
	}
	if(!validate_parameters(dist)){ goto cleanup; }

	return dist;

cleanup:
	free_Distribution(dist);
	return NULL;
}


int nparameter_distribution(const Distribution dist){
	return dist->np;
}

bool validate_parameters(const Distribution dist){
	int np = nparameter_distribution(dist);
	switch(dist->key){
		case 0:   break;
		case 'W': if(np!=2 || NULL==dist->param){ return false; }
			  if(dist->param[0]<0.){ return false; }  //shape
		          if(dist->param[1]<0.){ return false; }  //scale
			  break;
		case 'L': if(np!=2 || NULL==dist->param){ return false; }
			  if(dist->param[1]<0.){ return false; }  //scale
			  break;
		case 'N': if(np!=2 || NULL==dist->param){ return false; }
			  if(dist->param[1]<0.){ return false; }  //sd
			  break;
		case 'M': if(NULL==dist->info){ return false; }
			  NormMixParam normmix = (NormMixParam)dist->info;
			  if(np!=3*normmix->nmix-1){ return false; }
			  real_t sum=0.0;
			  for ( int i=0 ; i<normmix->nmix ; i++){
			  	if(normmix->sd[i]<0.){ return false; }  //sd
				if(normmix->prob[i]<0.|| normmix->prob[i]>1.){ return false; }  //prob
				sum += normmix->prob[i];
			  }
			  for ( int i=0 ; i<normmix->nmix ; i++){
				  normmix->prob[i] /= sum;
			  }
			  break;
		default: errx(EXIT_FAILURE,"Unrecognised distribution in %s",__func__);
	}
	return true;
}

real_t qdistribution(const real_t px, const Distribution dist, const bool tail, const bool logp){
	if(NULL==dist){ return NAN; }
	switch(dist->key){
		case 0:   return 0;
		case 'W': return qweibull(px,dist->param[0],dist->param[1],tail,logp);
		case 'L': return qlogistic(px,dist->param[0],dist->param[1],tail,logp);
	        case 'N': return qnorm(px,dist->param[0],dist->param[1],tail,logp);
		case 'M': return qmixnorm(px,(NormMixParam)dist->info,tail,logp);
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
