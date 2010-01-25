#include <string.h>
#include <tgmath.h>
#include "intensities.h"
#include "nuc.h"
#include "normal.h"

MODEL new_MODEL(const char * label, const uint32_t lane, const uint32_t tile, const real_t shape, const real_t scale, const MAT cov1, const MAT cov2){
    validate(shape>0,NULL);
    validate(scale>0,NULL);
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
    
    model->shape = shape;
    model->scale = scale;
    model->lane = lane;
    model->tile = tile;
    
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
    newmodel->lane      = model->lane;
    newmodel->tile      = model->tile;
    newmodel->shape     = model->shape;
    newmodel->scale     = model->scale;
    newmodel->paired    = model->paired;

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

void show_MODEL( FILE * fp, const MODEL model){
    validate(NULL!=fp,);
    validate(NULL!=model,);
    if(NULL!=model->label){fputs(model->label,fp);}
    fprintf(fp,"Parameters for %u cycle model\n",model->ncycle);
    fprintf(fp,"Lane %u, tile %u\n",model->lane,model->tile);
    fprintf(fp,"Brightness:\tshape=%f\tscale=%f\n",model->shape,model->scale);
    fputs("Covariance matrix of errors:\n",fp);
    show_MAT(fp,model->cov1,5,5);

    if(NULL!=model->cov2){
        fputs("Model is paired\n",fp);
        show_MAT(fp,model->cov2,5,5);
    }
}

MODEL new_MODEL_from_fp( FILE * fp){
    validate(NULL!=fp,NULL);
    char * label = NULL;
    char c = fgetc(fp); ungetc(c,fp);
    if ( '#'==c){
        size_t len = 0;
        char * ln = fgetln(fp,&len);
        label = calloc(len+1,sizeof(char));
        memcpy(label,ln,len*sizeof(char));
    }
    unsigned int lane,tile;
    fscanf(fp,"%u %u",&lane,&tile);
    real_t shape=0.0, scale=0.0;
    fscanf(fp,real_format_str real_format_str,&shape,&scale);
    MAT cov1 = new_MAT_from_fp(fp);
    MAT cov2 = new_MAT_from_fp(fp);
    MODEL model = new_MODEL(label,lane,tile,shape,scale,cov1,cov2);
    
    if(NULL!=cov2){free_MAT(cov2);};
    free_MAT(cov1);
    safe_free(label);
    return model;
}

MODEL trim_MODEL(const uint32_t ncycle, const MODEL model){
	MODEL newmod = NULL;
	MAT trimmedVar1=NULL, trimmedVar2=NULL;
	validate(NULL!=model,NULL);
	if(model->ncycle<ncycle){ return NULL; }
	trimmedVar1 = copy_MAT(trim_MAT(model->cov1,NBASE*ncycle,NBASE*ncycle,true));
	if(NULL==trimmedVar1){ goto cleanup; }
	trimmedVar2 = copy_MAT(trim_MAT(model->cov2,NBASE*ncycle,NBASE*ncycle,true));
        if(NULL!=model->cov2 && NULL==trimmedVar2){ goto cleanup; }
	newmod = new_MODEL(model->label,model->lane,model->tile,model->shape,model->scale,trimmedVar1,trimmedVar2);
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


MAT generate_pure_intensities ( const real_t lambda, const NUC * seq, const uint32_t ncycle, const MAT chol, MAT ints){
    validate(NULL!=seq,NULL);
    validate(NULL!=chol,NULL);
    if(NULL==ints){
        ints = new_MAT(NBASE*ncycle,1);
        validate(NULL!=ints,NULL);
    }
    
    rmultinorm(NULL,chol,NBASE*ncycle,ints);
    reshape_MAT(ints,NBASE);
    for ( uint32_t i=0 ; i<ncycle ; i++){
        // Assume that ambiguity in the sequence is a no call
        if ( NUC_AMBIG != seq[i]){
            ints->x[NBASE*i+seq[i]] += lambda;
        }
    }
    
    return ints;
}

real_t dchisq4( const real_t x){
    return exp(-0.5*x);
}
    
real_t lss ( const MAT x, const MAT invVchol ){
    validate(NULL!=x,NAN);
    validate(NULL!=invVchol,NAN);
    validate(invVchol->nrow==invVchol->ncol,NAN);
    validate(x->nrow== invVchol->ncol,NAN);

    // Form x^t invVchol invVchol^t x, noting that invVchol is lower diagonal
    const uint32_t n=x->nrow;    
    real_t tot = 0.;
    for ( uint32_t i=0 ; i<n ; i++){
        real_t s = 0.;
        for ( uint32_t j=0 ; j<=i ; j++){
            s += x->x[j] * invVchol->x[i*n+j];
        }
        tot += s*s;
    }
    return tot;
}
    
    
MAT likelihood_cycle_intensities ( const real_t lambda, const MAT ints, const MAT * invchol, MAT like){
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
            like->x[i*NBASE+j] = dchisq4(lss(tmp,invchol[i]));
            //like->x[i*NBASE+j] = dmultinorm(tmp,NULL,invchol[i],NBASE,true);
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

#ifdef TESTINT
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