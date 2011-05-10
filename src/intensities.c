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

#include <string.h>
#include <tgmath.h>
#include <stdint.h>
#include <inttypes.h>
#include "intensities.h"
#include "random.h"
#include "nuc.h"
#include "normal.h"
#include "lambda_distribution.h"

#define MODEL_FILE_VERSION 5

MODEL new_MODEL(const char * label, const char dist1, const real_t *param1, const char dist2, const real_t * param2, const MAT cov1, const MAT cov2){
    validate(NULL!=cov1,NULL);
    MODEL model = calloc(1,sizeof(*model));
    model->paired = false;
    
    if(!is_square(cov1)){ goto cleanup;} // Matrix not symmetric
    if( (cov1->nrow % NBASE) != 0 ){ goto cleanup;} // Matrix dimensions not compatible with intensity file
    const uint32_t ncycle = cov1->nrow/NBASE;
    model->ncycle = ncycle;
    model->orig_ncycle = ncycle;
    if ( NULL!=label){
        model->label = calloc(strlen(label)+1,sizeof(char));
        strcpy(model->label,label);
    }
    
    model->cov1 = copy_MAT(cov1);
    if(NULL==model->cov1){ goto cleanup; }
    model->chol1 = copy_MAT(cov1);
    if(NULL==model->chol1){ goto cleanup; }
    cholesky(model->chol1);
    // Split covariance into cycle
    model->invchol1 = block_diagonal_MAT(model->cov1,NBASE);
    if(NULL==model->invchol1){ goto cleanup; }
    for ( uint32_t i=0 ; i<ncycle ; i++){
        cholesky(model->invchol1[i]);
        invert_cholesky(model->invchol1[i]);
    }
    
    if( NULL!=cov2 ){
        if( (cov2->ncol!=cov1->ncol) || (cov2->nrow!=cov1->nrow) ){ goto cleanup;} // Matrices not 
        model->paired = true;
        model->cov2 = copy_MAT(cov2);
        if(NULL==model->cov2){ goto cleanup; }
        
        model->chol2 = copy_MAT(cov2);
        if(NULL==model->chol2){ goto cleanup; }
        cholesky(model->chol2);
        
        model->invchol2 = block_diagonal_MAT(model->cov2,NBASE);
        if(NULL==model->invchol2){ goto cleanup; }
        for ( uint32_t i=0 ; i<ncycle ; i++){
            cholesky(model->invchol2[i]);
            invert_cholesky(model->invchol2[i]);
        }
    }
    
    model->dist1 = dist1;
    int np = nparameter_distribution(dist1);
    model->param1 = calloc(np,sizeof(real_t));
    if(NULL==model->param1){ goto cleanup; }
    memcpy(model->param1,param1,np*sizeof(real_t));

    model->dist2 = dist2;
    np = nparameter_distribution(dist2);
    model->param2 = calloc(np,sizeof(real_t));
    if(NULL==model->param2){ goto cleanup; }
    memcpy(model->param2,param2,np*sizeof(real_t));
    
    return model;
    
cleanup:
    free_MODEL(model);
    return NULL;
}

void free_MODEL( MODEL model){
    validate(NULL!=model,);
    free_MAT(model->cov1);
    free_MAT(model->chol1);
    if(NULL!=model->invchol1){
        for( uint32_t i=0 ; i<model->orig_ncycle ; i++){
            free_MAT(model->invchol1[i]);
        }
        safe_free(model->invchol1);
    }
    free(model->param1);
    free(model->param2);
    free_MAT(model->cov2);
    free_MAT(model->chol2);
    if(NULL!=model->invchol2){
        for( uint32_t i=0 ; i<model->orig_ncycle ; i++){
            free_MAT(model->invchol2[i]);
        }
        safe_free(model->invchol2);
    }
    safe_free(model->label);
    safe_free(model);
}

MODEL copy_MODEL( const MODEL model){
    validate(NULL!=model,NULL)
    MODEL newmodel      = calloc(1,sizeof(*newmodel));
    newmodel->ncycle    = model->ncycle;
    newmodel->orig_ncycle = model->orig_ncycle;
    newmodel->dist1     = model->dist1;
    newmodel->dist2     = model->dist2;
    newmodel->paired    = model->paired;

    int np = nparameter_distribution(model->dist1);
    newmodel->param1 = calloc(np,sizeof(real_t));
    if(NULL==newmodel->param1){ goto cleanup; }
    memcpy(newmodel->param1,model->param1,np*sizeof(real_t));
    np = nparameter_distribution(newmodel->dist2);
    newmodel->param2 = calloc(np,sizeof(real_t));
    if(NULL==newmodel->param2){ goto cleanup; }
    memcpy(newmodel->param2,model->param2,np*sizeof(real_t));

    newmodel->cov1       = copy_MAT(model->cov1);
    if(NULL==newmodel->cov1){ goto cleanup; }
    newmodel->chol1      = copy_MAT(model->chol1);
    if(NULL==newmodel->chol1){ goto cleanup; }
    newmodel->invchol1    = calloc(model->orig_ncycle,sizeof(*newmodel->invchol1));
    if(NULL==newmodel->invchol1){ goto cleanup; }
    for ( uint32_t i=0 ; i<model->orig_ncycle ; i++){
        newmodel->invchol1[i]   = copy_MAT(model->invchol1[i]);
    }
    if(NULL==newmodel->invchol1){ goto cleanup; }
    
    newmodel->cov2       = copy_MAT(model->cov2);
    if(NULL==newmodel->cov2){ goto cleanup; }
    newmodel->chol2      = copy_MAT(model->chol2);
    if(NULL==newmodel->chol2){ goto cleanup; }
    newmodel->invchol2    = calloc(model->orig_ncycle,sizeof(*newmodel->invchol2));
    if(NULL==newmodel->invchol2){ goto cleanup; }
    for ( uint32_t i=0 ; i<model->orig_ncycle ; i++){
        newmodel->invchol2[i]   = copy_MAT(model->invchol2[i]);
    }
    if(NULL==newmodel->invchol2){ goto cleanup; }
    
    if(NULL!=model->label){
        newmodel->label = calloc(1+strlen(model->label),sizeof(char));
        strcpy(newmodel->label,model->label);
    }
    
    return newmodel;
    
cleanup:
    free_MODEL(newmodel);
    return NULL;
}

void show_MODEL( FILE * fp, MODEL model){
    validate(NULL!=fp,);
    validate(NULL!=model,);
    if(NULL!=model->label){fputs(model->label,fp);}
    fprintf(fp,"Parameters for %u cycle model\n",model->ncycle);
    fprintf(fp,"Brightness distribution:\n\tEnd1:%c",model->dist1);
    for ( int i=0 ; i<nparameter_distribution(model->dist1) ; i++){
        fprintf(fp,"\t%f",model->param1[i]);
    }
    if(NULL!=model->cov2){
	    fprintf(fp,"\tEnd2:%c",model->dist2);
	    for ( int i=0 ; i<nparameter_distribution(model->dist2) ; i++){
		    fprintf(fp,"\t%f",model->param2[i]);
	    }
	    fputc('\n',fp);
    } else {
	    fputc('\n',fp);
    }

    fputs("Covariance matrix of errors:\n",fp);
    show_MAT(fp,model->cov1,5,5);

    if(NULL!=model->cov2){
        fputs("Model is paired\n",fp);
        show_MAT(fp,model->cov2,5,5);
    }
}

MODEL new_MODEL_from_fp( FILE * fp){
    validate(NULL!=fp,NULL);

    char c = EOF;
    char * label = NULL;
    size_t labellen = 0;
    // Parse comments and forget
    while ( (c=fgetc(fp)) == '#' ){
        size_t len = 0;
        char * ln = NULL;
        #ifdef  _GNU_SOURCE
        getline(&ln,&len,fp);
        #else
        ln = fgetln(fp,&len);
        #endif

        label = reallocf(label,(1+labellen+len)*sizeof(char));
        if(NULL==label){ return NULL;} // Should really call cleanup here
        memcpy(label+labellen,ln,len*sizeof(char));
        label[labellen+len] = '\0';
        labellen += len;
        #ifdef _GNU_SOURCE
        free(ln);
        #endif
    }
    ungetc(c,fp);
    // Check version number
    int version=-1;
    fscanf(fp,"%*s %d",&version);
    if(version!=MODEL_FILE_VERSION){
	    errx(EXIT_FAILURE,"Version of model file is wrong. Expecting %d but got %d.",MODEL_FILE_VERSION,version);
    }

    uint32_t ncycle=0;
    char dist=0;
    int ret = fscanf(fp, "%"SCNu32 " %c",&ncycle,&dist);
    int np = nparameter_distribution(dist);
    real_t * param = calloc(np,sizeof(real_t));
    for ( int i=0 ; i<np ; i++){
    	int ret = fscanf(fp,real_format_str,&param[i]);
	if(ret!=1){
		errx(EXIT_FAILURE,"Problem parsing distribution parameter %d from model file.",i+1);
	}
    }
    if(!validate_parameters(param,dist)){
	    errx(EXIT_FAILURE,"Distribution parameters invalid");
    }

    MAT cov1 = new_MAT_from_fp(fp,ncycle*NBASE,ncycle*NBASE);


    uint32_t ncycle2=0;
    char dist2 = '\0';
    real_t * param2 = NULL;
    ret = fscanf(fp, "%"SCNu32 " %c",&ncycle2,&dist2);
    if( ret>0){
	    if(ret!=2 ){ errx(EXIT_FAILURE,"Problems reading second end of runfile: read %d elts",ret); }
    	    if(ncycle!=ncycle2){ errx(EXIT_FAILURE,"Reading paired-end runfile but ends have differing numbers of cycles (%"SCNu32 " and %"SCNu32 ")",ncycle,ncycle2); }
	    int np2 = nparameter_distribution(dist2);
	    param2 = calloc(np2,sizeof(real_t));
	    for ( int i=0 ; i<np2 ; i++){
		    int ret = fscanf(fp,real_format_str,&param2[i]);
		    if(ret!=1){
			    errx(EXIT_FAILURE,"Problem parsing distribution parameter %d from model file.",i+1);
		    }
	    }
	    if(!validate_parameters(param2,dist)){
		    errx(EXIT_FAILURE,"Distribution parameters invalid");
	    }
    }
    MAT cov2 = new_MAT_from_fp(fp,ncycle*NBASE,ncycle*NBASE);
    MODEL model = new_MODEL(label,dist,param,dist2,param2,cov1,cov2);
    
    if(NULL!=cov2){free_MAT(cov2);};
    free(param);
    free(param2);
    free_MAT(cov1);
    safe_free(label);
    return model;
}

MODEL trim_MODEL(const uint32_t ncycle, real_t final_factor[4], const MODEL model){
	MODEL newmod = NULL;
	MAT trimmedVar1=NULL, trimmedVar2=NULL;
	validate(NULL!=model,NULL);
	if(model->ncycle<ncycle){ return NULL; }
	trimmedVar1 = copy_MAT(trim_MAT(model->cov1,NBASE*ncycle,NBASE*ncycle,true));
	if(NULL==trimmedVar1){ goto cleanup; }
	trimmedVar2 = copy_MAT(trim_MAT(model->cov2,NBASE*ncycle,NBASE*ncycle,true));
        if(NULL!=model->cov2 && NULL==trimmedVar2){ goto cleanup; }

	// Alter covariance matrices
	if(final_factor[0]==-1.0){
		if(ncycle==model->ncycle){
			// If number of cycles match, don't scale
			final_factor[0] = final_factor[1] = final_factor[2] = final_factor[3] = 1.0;
		} else {
			// Learn scaling
			const size_t off1 = model->ncycle * 4 - 4;
			const size_t off2 = model->ncycle * 4 - 8;
			for ( int i=0 ; i<4 ; i++){
				final_factor[i] = model->cov1->x[(off1+i)*4*model->ncycle + off1+i] / model->cov1->x[(off2+i)*4*model->ncycle + off2+i];
				if(NULL!=model->cov2){
					final_factor[i] += model->cov2->x[(off1+i)*4*model->ncycle + off1+i] / model->cov2->x[(off2+i)*4*model->ncycle + off2+i];
					final_factor[i] /= 2;
				}
			}
		}
	}
	// Scale matrices
	fprintf(stderr,"Scaling variance of final cycle by %f %f %f %f\n",final_factor[0] , final_factor[1] , final_factor[2] , final_factor[3] );
	const size_t off = (ncycle * 4 - 4);
	for ( int i=0 ; i<4 ; i++){
		final_factor[i] = sqrt(final_factor[i]);
		for ( int j=0 ; j<ncycle*4 ; j++){
			trimmedVar1->x[(off+i)*4*ncycle+j] *= final_factor[i];
			trimmedVar1->x[j*ncycle*4+off+i] *= final_factor[i];
			if(NULL!=trimmedVar2){
				trimmedVar2->x[(off+i)*4*ncycle+j] *= final_factor[i];
				trimmedVar2->x[off+i+4*ncycle*j] *= final_factor[i];
			}
		}
	}

	newmod = new_MODEL(model->label,model->dist1,model->param1,model->dist2,model->param2,trimmedVar1,trimmedVar2);
	if(NULL==newmod){goto cleanup;}

	free_MAT(trimmedVar2);
	free_MAT(trimmedVar1);
	return newmod;

cleanup:
	free_MODEL(newmod);
	free_MAT(trimmedVar2);
	free_MAT(trimmedVar1);
	return NULL;
}

MODEL new_MODEL_from_file( const CSTRING filename ){
    FILE * fp = fopen(filename,"r");
    validate(NULL!=fp,NULL);
    MODEL model = new_MODEL_from_fp( fp);
    fclose(fp);
    return model;
}


MAT generate_pure_intensities ( 
    const real_t sdfact, const real_t lambda, const ARRAY(NUC) seq, 
    const ARRAY(NUC) adapter, const uint32_t ncycle, const MAT chol, 
    const real_t dustProb, const MAT invM, const MAT invP, const MAT N, MAT ints){
    validate(NULL!=seq.elt,NULL);
    validate(NULL!=chol,NULL);
    if(NULL==ints){
        ints = new_MAT(NBASE*ncycle,1);
        validate(NULL!=ints,NULL);
    }
    
    rmultinorm(NULL,chol,NBASE*ncycle,ints);
    reshape_MAT(ints,NBASE);
    if(1.0!=sdfact){scale_MAT(ints,sdfact);}
    for ( uint32_t i=0 ; i<ncycle ; i++){
        // Assume that ambiguity in the sequence is a no call
        if ( i<seq.nelt && NUC_AMBIG != seq.elt[i]){
            ints->x[NBASE*i+seq.elt[i]] += lambda;
        } else {
            uint32_t j = i-seq.nelt;
            if(j>=0 && j<adapter.nelt && NUC_AMBIG != adapter.elt[j]){
                ints->x[NBASE*i+adapter.elt[j]] += lambda;
            }
        }
    }

    // Add dust if required.
    if(0.0!=dustProb){
        real_t u = runif();
	if(u<dustProb){
	    const int cy = (int)(ncycle * u / dustProb);
	    const real_t dustval = lambda*10.0 - N->x[cy*NBASE+1];
	    for(int i=0 ; i<ncycle ; i++){
	        for ( int j=0 ; j<NBASE ; j++){
		    ints->x[i*NBASE+j] += dustval * invM->x[4+j] * invP->x[i*ncycle+cy];
		}
	    }
	    fprintf(stderr,"Dust at cycle %d\n",cy);
	}
    }
    
    return ints;
}

real_t dchisq4( const real_t x, const bool logb){
    return (false==logb)? exp(-0.5*x) : -0.5*x;
}
    
real_t lss ( const MAT x, const MAT invVchol ){
    validate(NULL!=x,NAN);
    validate(NULL!=invVchol,NAN);
    validate(invVchol->nrow==invVchol->ncol,NAN);
    validate(x->nrow==invVchol->ncol,NAN);

    // Form x^t invVchol invVchol^t x, noting that invVchol is lower triangular 
    // But lower triangle is stored in both the upper and lower triangle
    const uint32_t n=x->nrow;    
    real_t tot = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        real_t s = 0.;
        for ( uint32_t j=i ; j<n ; j++){
            s += x->x[j] * invVchol->x[i*n+j];
        }
        tot += s*s;
    }
    return tot;
}
    
    
MAT likelihood_cycle_intensities ( const real_t sdfact, real_t mu, const real_t lambda, const MAT ints, const MAT * invchol, MAT like){
    validate(NULL!=ints,NULL);
    validate(NULL!=invchol,NULL);
    const uint32_t ncycle = ints->ncol;
    
    if(NULL==like){
        like = new_MAT(NBASE,ncycle);
        validate(NULL!=like,NULL);
    }
    
    MAT tmp = new_MAT(NBASE,1);
    for ( int i=0 ; i<ncycle ; i++){
        for ( int j=0 ; j<NBASE ; j++){
            memcpy(tmp->x,ints->x+i*NBASE,NBASE*sizeof(real_t));
            tmp->x[j] -= lambda;
            like->x[i*NBASE+j] = -dchisq4(lss(tmp,invchol[i])/(sdfact*sdfact),true);
        }
    }
    if(mu>0.0){
        const real_t log_mu = log(mu);
        for( uint32_t i=0 ; i<(ncycle*NBASE) ; i++){
            // Calculate -log( mu + exp(loglike) )
            // == log(mu) + log( 1 + exp(loglike-log(mu)) )
            like->x[i] = -log_mu - log1p( exp(-like->x[i]-log_mu) );
        }
    }
    free_MAT(tmp);
    
    return like;
}

void fprint_vector(FILE * fp, const CSTRING prefix, const CSTRING sep, const CSTRING suffix, const real_t * x, const uint32_t n){
    validate(NULL!=fp,);
    validate(NULL!=prefix,);
    validate(NULL!=sep,);
    validate(NULL!=suffix,);
    validate(NULL!=x,);

    fputs(prefix,fp);
    if ( n>0){
        fprintf(fp,"%e",x[0]);
        for ( uint32_t i=1 ; i<n ; i++){
            fprintf(fp,"%s%e",sep,x[i]);
        }
    }
    fputs(suffix,fp);
}

// Print intensities as tab separated quad.
void fprint_intensities(FILE * fp, const char * prefix, const MAT ints, const bool last){
    validate(NULL!=fp,);
    validate(NULL!=ints,);
    validate(NBASE==ints->nrow,);
    
    if(NULL!=prefix){
        fputs(prefix,stdout);
        fprint_vector(fp,"\t"," ","",ints->x,NBASE);
    } else {
        fprint_vector(fp,""," ","",ints->x,NBASE);
    }
    const uint32_t ncycle = ints->ncol;
    for ( uint32_t cy=1 ; cy<ncycle ; cy++){
        fprint_vector(fp,"\t"," ","",ints->x+cy*NBASE,NBASE);
    }
    if(last){fputc('\n',fp);}
}

ARRAY(NUC) call_by_maximum_likelihood(const MAT likelihood, ARRAY(NUC) calls){
    validate(NULL!=likelihood,null_ARRAY(NUC));
    validate(NBASE==likelihood->nrow,null_ARRAY(NUC));
    const uint32_t ncycle = likelihood->ncol;
    if ( NULL==calls.elt){
       calls = new_ARRAY(NUC)(ncycle);
       validate(NULL!=calls.elt,null_ARRAY(NUC));
    }

    // likelihoods are stored as -log-likelihood so
    // max likelihood <==> min -log-likelihood
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        calls.elt[cycle] = NUC_A;
        real_t lmin = likelihood->x[cycle*NBASE];
        for ( uint32_t base=1 ; base<NBASE ; base++){
            if(likelihood->x[cycle*NBASE+base]<lmin){
                calls.elt[cycle] = base;
                lmin = likelihood->x[cycle*NBASE+base];
            }
        }
    }

    return calls;
}

real_t fuzz_prob(real_t prob){
	real_t scale = 1.0;
	real_t loc = 0.0;
	real_t nscale = 5.0;
	real_t l = 2.0 * scale * atanh(2.0*prob-1.0) + loc;
	l += nscale * (runif()-0.5);
	return 0.5*(1.0+tanh((l-loc)/(2.0*scale)));
}

ARRAY(PHREDCHAR) quality_from_likelihood(const MAT likelihood, const ARRAY(NUC) calls, const real_t generr, const bool doIllumina, ARRAY(PHREDCHAR) quals){
    validate(NULL!=likelihood,null_ARRAY(PHREDCHAR));
    validate(NULL!=calls.elt,null_ARRAY(PHREDCHAR));
    validate(NBASE==likelihood->nrow,null_ARRAY(PHREDCHAR));
    const uint32_t ncycle = likelihood->ncol;
    if ( NULL==quals.elt){
       quals = new_ARRAY(PHREDCHAR)(ncycle);
       validate(NULL!=quals.elt,null_ARRAY(PHREDCHAR));
    }
    
    // likelihoods are stored as -log-likelihood so
    // max likelihood <==> min -log-likelihood
    for ( uint32_t cycle=0 ; cycle<ncycle ; cycle++){
        NUC mb = calls.elt[cycle];
        real_t tot = 0., ml=likelihood->x[cycle*NBASE+mb];
        for( uint32_t b=0 ; b<NBASE ; b++){
            tot += exp(ml-likelihood->x[cycle*NBASE+b]);
        }
	real_t prob = (1.0-generr)/tot;
	//prob = fuzz_prob(prob);
        quals.elt[cycle] = phredchar_from_prob(prob,doIllumina);
    }
    return quals;
}


real_t purity ( real_t * ints4 ){
    real_t max1=fabs(ints4[0]);
    real_t max2=fabs(ints4[1]);
    if(max1<max2){ SWAP(max1,max2); }
    for ( int i=2 ; i<NBASE ; i++){
        real_t a = fabs(ints4[i]);
        if(max1<a){
            max2 = max1;
            max1 = a;
        } else if(max2<a){
            max2 = a;
        }
    }
    return max1/(max1+max2);
}

uint32_t number_inpure_cycles( const MAT intensities, const real_t threshold, const uint32_t ncycle){
    validate(NULL!=intensities,0);
    const real_t mcycle = (ncycle<intensities->ncol) ? ncycle : intensities->ncol;
    uint32_t count = 0;
    for ( uint32_t i=0 ; i<mcycle ; i++){
        real_t p = purity(intensities->x+i*NBASE);
        if(p<threshold){ count++;}
    }
    return count;
}

/**
 * Produce raw intensities from simulated processed intensities.
 * Based on process_intensities from AYB
 */
MAT unprocess_intensities(const MAT intensities, const MAT M_t, const MAT P_t, const MAT N, MAT p){
    validate(NULL!=intensities,NULL);
    validate(NULL!=M_t,NULL);
    validate(NULL!=P_t,NULL);
    validate(NULL!=N,NULL);
    
    const uint_fast32_t ncycle = P_t->nrow;
    if(NULL==p){
        p = new_MAT(NBASE,ncycle);
        validate(NULL!=p,NULL);
    }
    bzero(p->x,p->nrow * p->ncol * sizeof(real_t));

    for( uint_fast32_t icol=0 ; icol<ncycle ; icol++){    // Columns of Intensity
        real_t dp[NBASE] = {0,0,0,0};
        for( uint_fast32_t chan=0 ; chan<NBASE ; chan++){ // Channels
            for ( uint_fast32_t base=0 ; base<NBASE ; base++){  // Bases (cols of M, rows of Minv_t)
                dp[chan] += M_t->x[chan*NBASE+base] * intensities->x[icol*NBASE+base];
            }
        }
        for ( uint_fast32_t pcol=0 ; pcol<ncycle ; pcol++){ // Columns of p
            const real_t tmp = P_t->x[icol*ncycle+pcol];
            for( uint_fast32_t chan=0 ; chan<NBASE ; chan++){
                p->x[pcol*NBASE+chan] += tmp * dp[chan];
            }
        }
        for ( uint_fast32_t chan=0 ; chan<NBASE ; chan++){
            p->x[icol*NBASE+chan] += N->x[icol*NBASE+chan];
        }
    }
    return p;
}


#ifdef TEST
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "random.h"


int main(int argc, char * argv[]){
    if(argc!=6){
        fputs("Usage: test lambda seq varfile n seed\n",stderr);
        return EXIT_FAILURE;
    }

    real_t lambda=0;
    sscanf(argv[1], real_format_str, &lambda);

    const unsigned int seqlen = strlen(argv[2]);
    NUC * seq = calloc(seqlen,sizeof(NUC));
    for( int i=0 ; i<seqlen ; i++){
        seq[i] = nuc_from_char(argv[2][i]);
    }

/*    MAT var = read_MAT_from_file(argv[3]);
    if(var->nrow!=var->ncol){
        fprintf(stderr,"Invalid variance file (%d row, %d col)\n",var->nrow,var->ncol);
        return EXIT_FAILURE;
    }
    if(var->nrow!=4*seqlen){
        fprintf(stderr,"Sequence length is %u but matrix has %d rows (should be %d).\n",seqlen,var->nrow,4*var->nrow);
        return EXIT_FAILURE;
    }
*/
    MAT chol = cholesky(identity_MAT(4*seqlen));
    MAT invchol = invertcholesky(copy_MAT(chol));
    
    unsigned int n = 0;
    sscanf(argv[4],"%u",&n);

    long unsigned int seed = 0;
    sscanf(argv[5],"%lu",&seed);
    init_gen_rand(seed);
    
    fprintf(stdout,"Generating %u intensities for %u cycles.\n",n,seqlen);
    MAT ints = NULL;
    MAT like = NULL;
    for ( uint32_t i=0 ; i<n ; i++){
        fprintf(stdout,"* Set %d\n",i+1);
        ints = generate_pure_intensities( lambda, seq, seqlen, chol, ints);
        fputs("Intensities\n",stdout);
        show_MAT(stdout,ints,5,6);
        like = likelihood_cycle_intensities( lambda,ints,invchol,like);
        fputs("Log-likelihoods\n",stdout);
        show_MAT(stdout,like,5,6);
    }
    
    return EXIT_SUCCESS;
}
#endif
