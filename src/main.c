#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "sequence.h"
#include "random.h"
#include "utility.h"
#include "intensities.h"
#include "weibull.h"


void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"simNGS\"\n"
"Simulate likelihoods for Illumina data from fasta format files\n"
"\n"
"Usage:\n"
"\tsimNGS [-c ncycle] [-l shape,scale] [-p on|off] [-r mu]\n"
"\t       [-s seed] runfile\n"
"\tsimNGS --help\n"
"simNGS reads from stdin and writes to stdout, messages and progess\n"
"indicators are written to stderr.\n"
,fp);
}

void fprint_help( FILE * fp){
    validate(NULL!=fp,);
    fputs(
/*
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/
"\n"
"-c, --cycles ncycle [default: as runfile]\n"
"\tNumber of cycles to do, up to maximum allowed for runfile.\n"
"\n"
"-l, --lambda shape,scale [default: as runfile]\n"
"\tShape and scale of cluster brightness distribution (lambda).\n"
"Currently a Weibull distribution is used.\n"
"\n"
"-p, --paired on|off [default: as runfile]\n"
"\tTreat run as paired-end or not. For single-ended runs treated as\n"
"paired, the covariance matrix is duplicated to make two uncorrelated pairs.\n"
"For paired-end runs treated as single, the second end is ignored.\n"
"\n"
"-r, --robust mu [default: 0]\n"
"\tCalculate robustified likelihood, equivalent to adding mu to every\n"
"likelihood.\n"
"\n"
"-s, --seed seed [default: clock]\n"
"\tSet seed from random number generator.\n"
, fp);
}

static struct option longopts[] = {
    { "cycles",     required_argument, NULL, 'c' },
    { "lambda",     required_argument, NULL, 'l' },
    { "paired",     required_argument, NULL, 'p' },
    { "robust",     required_argument, NULL, 'r' },
    { "seed",       required_argument, NULL, 's' },
    { "help",       no_argument,       NULL, 'h' }
};

bool parse_bool( const CSTRING str){
    validate(NULL!=str,false);
    if ( toupper(str[0]) == 'T' ){ return true; }
    if ( toupper(str[0]) == 'F' ){ return false; }
    if ( strcmp(str,"on") == 0 ){ return true; }
    if ( strcmp(str,"off") == 0 ){ return false;}
    int val=0;
    int ret = sscanf(str,"%d",&val);
    if(ret==1 && val==1){ return true;}
    // Combining remaining and failure cases below
    return false;
}

const CSTRING boolstr[] = { "false", "true" };

real_t parse_real( const CSTRING str){
    validate(NULL!=str,NAN);
    real_t x = NAN;
    sscanf(str,real_format_str,&x);
    return x;
}

unsigned int parse_uint( const CSTRING str){
    validate(NULL!=str,0);
    unsigned int n=0;
    sscanf(str,"%u",&n);
    return n;
}

typedef struct {
    unsigned int ncycle;
    real_t shape,scale;
    bool paired,paired_set;
    real_t mu;
    uint32_t seed;
} * SIMOPT;

SIMOPT new_SIMOPT(void){
    SIMOPT opt = calloc(1,sizeof(*opt));
    validate(NULL!=opt,NULL);
    // Set defaults
    opt->ncycle = 0;
    opt->shape = 0.0;
    opt->scale = 0.0;
    opt->paired_set = false;
    opt->mu = 0.0;
    opt->seed = 0;
    
    return opt;
}

void free_SIMOPT(SIMOPT opt){
    validate(NULL!=opt,);
    safe_free(opt);
}

SIMOPT copy_SIMOPT ( const SIMOPT simopt){
    validate(NULL!=simopt,NULL);
    SIMOPT newopt = calloc(1,sizeof(*newopt));
    validate(NULL!=newopt,NULL);
    memcpy(newopt,simopt,sizeof(*simopt));
    return newopt;
}

void show_SIMOPT (FILE * fp, const SIMOPT simopt){
    validate(NULL!=fp,);
    validate(NULL!=simopt,);
    fputs("\tOptions:\n",fp);
    fprintf( fp,"ncycle\t%u\n",simopt->ncycle);
    fprintf( fp,"paired\t%s\n",boolstr[simopt->paired]);
    fprintf( fp,"mu\t%f\n",simopt->mu);
    fprintf( fp,"shape\t%f\n",simopt->shape);
    fprintf( fp,"scale\t%f\n",simopt->scale);
    fprintf( fp,"seed\t%u\n",simopt->seed);
}

SIMOPT parse_arguments( const int argc, char * const argv[] ){
    int ch;
    SIMOPT simopt = new_SIMOPT();
    validate(NULL!=simopt,NULL);
    
    while ((ch = getopt_long(argc, argv, "c:l:p:r:s:h", longopts, NULL)) != -1){
        switch(ch){
        case 'c':   sscanf(optarg,"%u",&simopt->ncycle);
                    break;
        case 'l':   sscanf(optarg,real_format_str "," real_format_str ,&simopt->shape,&simopt->scale);
                    break;
        case 'p':   simopt->paired = parse_bool(optarg);
                    simopt->paired_set = true;
                    break;
        case 'r':   simopt->mu = parse_real(optarg);
                    break;
        case 's':   simopt->seed = parse_uint(optarg);
                    break;
        case 'h':
            fprint_help(stderr);
            exit(EXIT_SUCCESS);
        default:
            fprint_usage(stderr);
            exit(EXIT_FAILURE);
        }
    }

    return simopt;
}

int main( int argc, char * argv[] ){
    SIMOPT simopt = parse_arguments(argc,argv);

    argc -= optind;
    argv += optind;
    if(0==argc){
        fputs("Expecting runfile on commandline but none found.\n",stderr);
        fprint_usage(stderr);
        return EXIT_FAILURE;
    }

    // Load up model
    MODEL model = new_MODEL_from_file(argv[0]);
    // Resolve options and model
    if(simopt->shape!=0){ model->shape = simopt->shape;}
    simopt->shape = model->shape;
    if(simopt->scale!=0){ model->scale = simopt->scale;}
    simopt->scale = model->scale;
    
    if(simopt->paired_set && (simopt->paired!=model->paired)){
        if(simopt->paired==false){
            model->paired = false;
            free_MAT(model->cov2);
            model->cov2 = NULL;
        } else {
            model->paired = true;
            model->cov2 = copy_MAT(model->cov1);
            model->chol2 = copy_MAT(model->chol1);
            model->invchol2 = calloc(model->ncycle,sizeof(*model->invchol2));
            for ( uint32_t i=0 ; i<model->ncycle ; i++){
                model->invchol2[i] = copy_MAT(model->invchol1[i]);
            }
        }
    }
    simopt->paired = model->paired;
    
    if(simopt->ncycle!=0){
        if(simopt->ncycle>model->ncycle){
            fprintf(stderr,"Asked for more cycles than runfile allows. Doing %u.\n",model->ncycle);
        } else {
            MODEL newmodel = trim_MODEL(simopt->ncycle,model);
            free_MODEL(model);
            model = newmodel; 
        }
    }
    simopt->ncycle = model->ncycle;
    
    
    // Initialise random number generator
    if ( simopt->seed==0 ){
        uint32_t seed = (uint32_t) time(NULL);
        fprintf(stderr,"Using seed %u\n",seed);
        simopt->seed = seed;
    }
    init_gen_rand( simopt->seed );
    //show_SIMOPT(stderr,simopt);
    //show_MODEL(stderr,model);


    // Scan through fasta file
    MAT intensities = NULL, intensities2 = NULL;
    MAT loglike = NULL, loglike2=NULL;
    SEQ seq = NULL;
    while ((seq=sequence_from_fasta(stdin))!=NULL){
        //show_SEQ(stderr,seq);
        if ( seq->length < model->ncycle ){
            fprintf(stderr,"Sequence %s shorter than number of cycles, skipping\n",seq->name);
            free_SEQ(seq);
            continue;
        }
        real_t lambda = rweibull(model->shape,model->scale);
        intensities = generate_pure_intensities(lambda,seq->seq,model->ncycle,model->chol1,intensities);
        loglike = likelihood_cycle_intensities(lambda,intensities,model->invchol1,loglike);
        uint32_t x = (uint32_t)( 1794 * runif());
        uint32_t y = (uint32_t)( 2048 * runif());
        fprintf(stdout,"%u\t%u\t%u\t%u",model->lane,model->tile,x,y);
        fprint_intensities(stdout,"",loglike,false);
        if ( model->paired ){
            NUC * rcseq = reverse_complement(seq->seq,seq->length);
            intensities2 = generate_pure_intensities(lambda,rcseq,model->ncycle,model->chol2,intensities2);
            loglike2 = likelihood_cycle_intensities(lambda,intensities2,model->invchol2,loglike2);
            fprint_intensities(stdout,"",loglike2,false);
            safe_free(rcseq);
        }
        fputc('\n',stdout);
        free_SEQ(seq);
    }
    return EXIT_SUCCESS;
}

