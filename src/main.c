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
#include "normal.h"


#define ILLUMINA_ADAPTER "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"

void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"simNGS\"\n"
"Simulate likelihoods for Illumina data from fasta format files\n"
"\n"
"Usage:\n"
"\tsimNGS [-a adapter] [-b shape:scale] [-c correlation] [-d]\n"
"\t       [-f nimpure:ncycle:threshold] [-i filename] [-l lane]\n"
"\t       [-m insertion:deletion:mutation] [-n ncycle] [-o output_format]\n"
"\t       [-p] [-q quantile] [-r mu] [-s seed] [-t tile] [-v factor ]\n"
"\t       runfile\n"
"\tsimNGS --help\n"
"\tsimNGS --licence\n"
"simNGS reads from stdin and writes to stdout. Messages and progess\n"
"indicators are written to stderr.\n"
"\n"
"Example:\n"
"\tcat sequences.fa | simNGS runfile > sequences.like\n"
,fp);
}

void fprint_licence (FILE * fp){
    validate(NULL!=fp,);
    fputs(
#define PROGNAME "simNGS"
"  simNGS software for simulating likelihoods for next-gen sequencing machines\n"
#include "copyright.inc"
    ,fp);
}

void fprint_help( FILE * fp){
    validate(NULL!=fp,);
    fputs(
/*
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/
"\n"
"-a, --adapter sequence [default: " ILLUMINA_ADAPTER " ]\n"
"\tSequence to pad reads with if they are shorter than the number of\n"
"cycles required, reflecting the adpater sequence used for sample preparaton.\n"
"A null adapter (i.e. pad with ambiguity characters) can be specified by -a \"\"\n" 
"\n" 
"-b, --brightness shape:scale [default: as runfile]\n"
"\tShape and scale of cluster brightnes distribution.\n"
"Currently a Weibull distribution is used.\n"
"\n"
"-c, --correlation [default: 1.0]\n"
"\tCorrelation between the cluster brightness of one end of a paired-end\n"
"run and the other. Default is complete correlation, the ends having equal\n"
"brightness. Correlation should belong to [-1,1].\n"
"\n"
"-d, --describe\n"
"\tPrint a description of the runfile and exit.\n"
"\n"
"-f, --filter nimpure:ncycle:threshold [default: no filtering]\n"
"\tUse purity filtering on generated intensities, allowing a maximum of\n"
"nimpure cyles in the first ncycles with a purity greater than threshold.\n"
"\n"
"-i, --intensities filename [default: none]\n"
"\tWrite the processed intensities generated to \"filename\".\n"
"\n"
"-l, --lane lane [default: as runfile]\n"
"\tSet lane number\n"
"\n"
"-m, --mutate insertion:deletion:mutation [default: no mutation]\n"
"\tSimple model of sequence mutation to reflect sample preparation errors.\n"
"Each sequence read in is transformed by a simple automata which inserts,"
"deletes or mutates bases with the specified probabilities.\n"
"\n"
"-n, --ncycles ncycles [default: as runfile]\n"
"\tNumber of cycles to do, up to maximum allowed for runfile.\n"
"\n"
"-o, --output format [default: likelihood]\n"
"\tFormat in which to output results. Either \"likelihood\" or \"fasta\"\n"
"\n"
"-p, --paired\n"
"\tTreat run as paired-end. For single-ended runs treated as\n"
"paired, the covariance matrix is duplicated to make two uncorrelated pairs.\n"
"For paired-end runs treated as single, the second end is ignored.\n"
"\n"
"-q, --quantile quantile [default: 0]\n"
"\tQuantile below which cluster brightness is discarded and redrawn from\n"
"distribution.\n"
"\n"
"-r, --robust mu [default: 0]\n"
"\tCalculate robustified likelihood, equivalent to adding mu to every\n"
"likelihood.\n"
"\n"
"-s, --seed seed [default: clock]\n"
"\tSet seed from random number generator.\n"
"\n"
"-t, --tile tile [default: as runfile\n"
"\tSet tile number.\n"
"\n"
"-v, --variance factor [default: 1.0]\n"
"\tFactor with which to scale variance matrix by.\n"
, fp);
}

static struct option longopts[] = {
    { "adapter",    required_argument, NULL, 'a'},
    { "brightness", required_argument, NULL, 'b' },
    { "correlation", required_argument, NULL, 'c' },
    { "describe",   no_argument,       NULL, 'd' },
    { "filter",     required_argument, NULL, 'f' },
    { "intensities", required_argument, NULL, 'i'},
    { "lane",       required_argument, NULL, 'l' },
    { "mutate",     required_argument, NULL, 'm' },
    { "ncycle",     required_argument, NULL, 'n' },
    { "output",     required_argument, NULL, 'o' },
    { "paired",     no_argument,       NULL, 'p' },
    { "quantile",   required_argument, NULL, 'q' },
    { "robust",     required_argument, NULL, 'r' },
    { "seed",       required_argument, NULL, 's' },
    { "tile",       required_argument, NULL, 't' },
    { "variance",   required_argument, NULL, 'v' },
    { "help",       no_argument,       NULL, 'h' },
    { "licence",    no_argument,       NULL, 0 },
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

enum outformat { OUTPUT_LIKE=0, OUTPUT_FASTA };
const char * output_format_str[] = { "likelihood", "fasta" };

typedef struct {
    unsigned int ncycle;
    real_t shape,scale,corr,threshold;
    bool paired,desc;
    real_t mu,sdfact;
    uint32_t seed;
    uint32_t tile,lane;
    real_t purity_threshold;
    uint32_t purity_cycles,purity_max;
    CSTRING intensity_fn;
    enum outformat format;
    bool mutate;
    real_t ins,del,mut;
    ARRAY(NUC) adapter;
} * SIMOPT;

SIMOPT new_SIMOPT(void){
    SIMOPT opt = calloc(1,sizeof(*opt));
    validate(NULL!=opt,NULL);
    // Set defaults
    opt->ncycle = 0;
    opt->shape = 0.0;
    opt->scale = 0.0;
    opt->corr = 1.0;
    opt->threshold = 0.0;
    opt->paired = false;
    opt->desc = false;
    opt->mu = 0.0;
    opt->seed = 0;
    opt->sdfact = 1.0;
    opt->tile = 0;
    opt->lane = 0;
    opt->purity_threshold = 0.;
    opt->purity_cycles = 0;
    opt->purity_max = 0;
    opt->intensity_fn = NULL;
    opt->format = OUTPUT_LIKE;
    opt->mutate = false;
    opt->ins=0.; opt->del=0.; opt->mut=0.;
    opt->adapter = nucs_from_string(ILLUMINA_ADAPTER);
    
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
    if(NULL!=simopt->intensity_fn){
        newopt->intensity_fn = copy_CSTRING(simopt->intensity_fn);
    }
    return newopt;
}

void show_SIMOPT (FILE * fp, const SIMOPT simopt){
    validate(NULL!=fp,);
    validate(NULL!=simopt,);
    fputs("\tOptions:\n",fp);
    fprintf( fp,"ncycle\t%u\n",simopt->ncycle);
    fprintf( fp,"paired\t%s\n",boolstr[simopt->paired]);
    fprintf( fp,"Brightness correlation\t%f\n",simopt->corr);
    fprintf( fp,"mu\t%f\n",simopt->mu);
    fprintf( fp,"shape\t%f\n",simopt->shape);
    fprintf( fp,"scale\t%f\n",simopt->scale);
    fprintf( fp,"threshold\t%f\n",simopt->threshold);
    fprintf( fp,"variance factor\t%f\n",simopt->sdfact*simopt->sdfact);
    fprintf( fp,"tile\t%u\tlane%u\n",simopt->tile,simopt->lane);
    fprintf( fp,"seed\t%u\n",simopt->seed);

    if(simopt->purity_cycles!=0){
       fprintf( fp,"Purity filtering: threshold %f. Maximum of %u inpure in %u cycles\n",simopt->purity_threshold, simopt->purity_max,simopt->purity_cycles);
    } else {
       fputs("No purity filtering.\n",fp);
    }

    if(simopt->ins!=0. && simopt->del!=0. && simopt->mut!=0.){
        fprintf(fp,"Input sequence will be mutated with ins %f, del %f, mut %f\n",simopt->ins,simopt->del,simopt->mut);
    }

    fputs("Reads will be padded with adapter sequence if necessary.\n",fp);
    show_ARRAY(NUC)(fp,simopt->adapter,"",0);
    fputc('\n',fp);

    fprintf(fp,"Writing %s to output.\n",output_format_str[simopt->format]);
    if(NULL!=simopt->intensity_fn){
        fprintf(fp,"Will write intensities to \"%s\"\n",simopt->intensity_fn);
    }
}

SIMOPT parse_arguments( const int argc, char * const argv[] ){
    int ch;
    SIMOPT simopt = new_SIMOPT();
    validate(NULL!=simopt,NULL);
    
    while ((ch = getopt_long(argc, argv, "a:b:c:df:i:l:m:n:o:pq:r:s:t:uv:h", longopts, NULL)) != -1){
        int ret;
        unsigned long int i=0,j=0;
        switch(ch){
        case 'a':   free_ARRAY(NUC)(simopt->adapter);
                    simopt->adapter = nucs_from_string(optarg);
                    break;
        case 'b':   ret = sscanf(optarg,real_format_str ":" real_format_str ,&simopt->shape,&simopt->scale);
                    if(ret!=2){ errx(EXIT_FAILURE,"Insufficient arguments for brightness.");}
                    if(simopt->shape<=0.){ errx(EXIT_FAILURE,"Brightness shape must be greater than zero."); }
                    if(simopt->scale<=0.){ errx(EXIT_FAILURE,"Brightness scale must be greater than zero."); }
                    break;
        case 'c':   sscanf(optarg,"%f",&simopt->corr);
                    if(simopt->corr<-1.0 || simopt->corr>1.0){errx(EXIT_FAILURE,"Correlation between end brightness should be in [-1,1]. Was given %f.",simopt->corr);}
                    break;
        case 'd':   simopt->desc = true;
                    break;
        case 'f':   ret = sscanf(optarg, "%lu:%lu:" real_format_str,&i,&j,&simopt->purity_threshold);
                    if(ret!=3){ errx(EXIT_FAILURE,"Insufficient arguments for filtering.");}
                    simopt->purity_max = i;
                    simopt->purity_cycles = j;
                    if( simopt->purity_threshold<0. || simopt->purity_threshold>1.0){
                        errx(EXIT_FAILURE,"Purity threshold is %f but should be between 0 and 1.",simopt->purity_threshold);
                    }
                    break;
        case 'i':   simopt->intensity_fn = copy_CSTRING(optarg);
                    break;
        case 'l':   simopt->lane = parse_uint(optarg);
                    if(simopt->lane==0){errx(EXIT_FAILURE,"Lane number must be greater than zero.");}
                    break;
        case 'm':   ret = sscanf(optarg, real_format_str ":" real_format_str ":" real_format_str,&simopt->ins,&simopt->del,&simopt->mut);
                    if( ret!=3 ){ errx(EXIT_FAILURE,"Insufficient arguments for mutation.");}
                    if(!isprob(simopt->ins) || !isprob(simopt->del) || !isprob(simopt->mut) ){
                        errx(EXIT_FAILURE,"Mutation parameters not probabilities. Given: ins %f, del %f, mut %f",simopt->ins,simopt->del,simopt->mut);
                    }
                    if(simopt->ins+simopt->del+simopt->mut>1.0){
                        errx(EXIT_FAILURE,"Mutation parameters sum to greater than one.");
                    }
                    simopt->mutate = true;
                    break;   
        case 'n':   sscanf(optarg,"%u",&simopt->ncycle);
                    if(simopt->ncycle==0){errx(EXIT_FAILURE,"Number of cycles to simulate must be greater than zero.");}
                    break;
        case 'o':   if( strcasecmp(optarg,output_format_str[OUTPUT_LIKE])==0 ){ simopt->format = OUTPUT_LIKE; }
                    else if ( strcasecmp(optarg,output_format_str[OUTPUT_FASTA])==0 ){ simopt->format = OUTPUT_FASTA; }
                    else {
                        errx(EXIT_FAILURE,"Unrecognised output option %s.",optarg);
                    }
                    break;
        case 'p':   simopt->paired = true;
                    break;
        case 'q':   simopt->threshold = parse_real(optarg);
                    if(!isprob(simopt->threshold) ){ 
                       errx(EXIT_FAILURE,"Threshold quantile to discard brightness must be a probability (got %e)\n",simopt->threshold);
                    }
                    break;
        case 'r':   simopt->mu = parse_real(optarg);
                    if(simopt->mu<0.0){errx(EXIT_FAILURE,"Robustness \"mu\" must be non-negative.");}
                    break;
        case 's':   simopt->seed = parse_uint(optarg);
                    break;
        case 't':   simopt->tile = parse_uint(optarg);
                    if(simopt->tile==0){errx(EXIT_FAILURE,"Tile number must be greater than zero.");}
                    break;
        case 'v':   simopt->sdfact = parse_real(optarg);
                    if(simopt->sdfact<0.0){errx(EXIT_FAILURE,"Variance scaling factor must be non-negative.");}
                    simopt->sdfact = sqrt(simopt->sdfact);
                    break;
        case 'h':
            fprint_usage(stderr);
            fprint_help(stderr);
            exit(EXIT_SUCCESS);
        case 0:
            fprint_licence(stderr);
            exit(EXIT_SUCCESS);
        default:
            fprint_usage(stderr);
            exit(EXIT_FAILURE);
        }
    }

    return simopt;
}

static inline real_t phred(const real_t p) __attribute__((const));
static inline real_t prop_upper( const real_t p, const uint32_t n) __attribute__((const));
static inline real_t prop_lower( const real_t p, const uint32_t n) __attribute__((const));

static inline real_t phred(const real_t p){
    return -10.0 * log10(p);
}

// Proportion confidence interval using Wilson's method
static inline real_t prop_upper( const real_t p, const uint32_t n){
    const real_t z = 1.959964;
    real_t desc = p*(1-p)/n + (z*z/(4.0*n*n));
    desc = p + z*z/(2.*n) + z * sqrt(desc);
    return desc / (1.0 + z*z/n);
}

static inline real_t prop_lower( const real_t p, const uint32_t n){
    const real_t z = 1.959964;
    real_t desc = p*(1-p)/n + (z*z/(4.0*n*n));
    desc = p + z*z/(2.*n) - z * sqrt(desc);
    return desc / (1.0 + z*z/n);
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
    if (NULL==model){
        errx(EXIT_FAILURE,"Failed to read runfile \"%s\"",argv[0]);
    }
    if( simopt->desc ){
        show_MODEL(stderr,model);
        return EXIT_SUCCESS;
    } 
    fprintf(stderr,"Description of runfile:\n%s",model->label);
    
    // Resolve options and model
    // Should factor out into separate routine
    if(simopt->shape!=0){ model->shape = simopt->shape;}
    simopt->shape = model->shape;
    if(simopt->scale!=0){ model->scale = simopt->scale;}
    simopt->scale = model->scale;
    
    if(simopt->paired!=model->paired){
        if(simopt->paired==false){
            fputs("Treating paired-end model as single-ended.\n",stderr);
            model->paired = false;
            free_MAT(model->cov2);
            model->cov2 = NULL;
        } else {
            fputs("Treating single-ended model as paired-end.\n",stderr);
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
    if(0!=simopt->lane){ model->lane = simopt->lane;}
    if(0!=simopt->tile){ model->tile = simopt->tile;}
    simopt->tile = model->tile;
    simopt->lane = model->lane;

    
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
    MAT intensities = NULL, loglike = NULL;
    ARRAY(NUC) calls = null_ARRAY(NUC);
    SEQ seq = NULL;
    FILE * fp =  stdin; //fopen("test/test100_small.fa","r");
    uint32_t seq_count=0, unfiltered_count=0;
    // Memory for error counting
    uint32_t * error = calloc(model->ncycle,sizeof(uint32_t));
    uint32_t * error2 = calloc(model->ncycle,sizeof(uint32_t));
    uint32_t errorhist[7] = {0,0,0,0,0,0,0};
    uint32_t errorhist2[7] = {0,0,0,0,0,0,0};
    
    FILE * fpout = (NULL!=simopt->intensity_fn) ? fopen(simopt->intensity_fn,"w") : NULL;
    if ( NULL==fpout && NULL!=simopt->intensity_fn){
        fprintf(stderr,"Failed to open \"%s\" for writing.\n",simopt->intensity_fn);
    }

    while ((seq=sequence_from_fasta(fp))!=NULL){
        //show_SEQ(stderr,seq);
        if(simopt->mutate){
            SEQ mut = mutate_SEQ(seq,simopt->ins,simopt->del,simopt->mut);
            free_SEQ(seq);
            seq = mut;
        }

        // Pick lambda using Gaussian Copula
        real_t lambda1=NAN,lambda2=NAN;
        {
            real_t px=0.0,py=0.0;
            do{
               const real_t corr = simopt->corr;
               // Two correlated Gaussians
               real_t x = rstdnorm();
               real_t y = corr*x + sqrt(1-corr*corr) * rstdnorm();
               // Convert to uniform deviates (the copula)
               px = pstdnorm(x,false,false);
               py = pstdnorm(y,false,false);
            } while(px<simopt->threshold || py<simopt->threshold);
            // Convert to Weibull via inversion formula
            lambda1 = qweibull(px,simopt->shape,simopt->scale,false,false);
            lambda2 = qweibull(py,simopt->shape,simopt->scale,false,false);
        }
            
        intensities = generate_pure_intensities(simopt->sdfact,lambda1,seq->seq,simopt->adapter,model->ncycle,model->chol1,intensities);
        loglike = likelihood_cycle_intensities(simopt->sdfact,simopt->mu,lambda1,intensities,model->invchol1,loglike);
        uint32_t x = (uint32_t)( 1794 * runif());
        uint32_t y = (uint32_t)( 2048 * runif());

        if(NULL!=fpout){
            fprintf(fpout,"%u\t%u\t%u\t%u",model->lane,model->tile,x,y);
            fprint_intensities(fpout,"",intensities,false);
        }

        switch(simopt->format){
           case OUTPUT_LIKE:
               fprintf(stdout,"%u\t%u\t%u\t%u",model->lane,model->tile,x,y);
               break;
           case OUTPUT_FASTA:
               fprintf(stdout,">%s\n",seq->name);
               break;
        }

        if ( number_inpure_cycles(intensities,simopt->purity_threshold,simopt->purity_cycles) <= simopt->purity_max){
            calls = call_by_maximum_likelihood(loglike,calls);
            uint32_t nerr = 0;
            for ( uint32_t i=0 ; i<model->ncycle ; i++){
                if(calls.elt[i] != seq->seq.elt[i]){ nerr++; error[i]++;}
            }
            errorhist[(nerr<6)?nerr:6]++;
            switch(simopt->format){
                case OUTPUT_LIKE:
                    fprint_intensities(stdout,"",loglike,false);
                    break;
                case OUTPUT_FASTA:
                    show_ARRAY(NUC)(stdout,calls,"",0);
                    break;
            }
            if ( model->paired ){
                ARRAY(NUC) rcseq = reverse_complement(seq->seq);
                intensities = generate_pure_intensities(simopt->sdfact,lambda2,rcseq,simopt->adapter,model->ncycle,model->chol2,intensities);
                loglike = likelihood_cycle_intensities(simopt->sdfact,simopt->mu,lambda2,intensities,model->invchol2,loglike);
                if(NULL!=fpout){ fprint_intensities(fpout,"",intensities,false); }
                calls = call_by_maximum_likelihood(loglike,calls);
                uint32_t nerr=0;
                for ( uint32_t i=0 ; i<model->ncycle ; i++){
                    if(calls.elt[i] != rcseq.elt[i]){ nerr++; error2[i]++;}
                }
                errorhist2[(nerr<6)?nerr:6]++;
                switch(simopt->format){
                    case OUTPUT_LIKE:
                        fprint_intensities(stdout,"",loglike,false);
                        break;
                    case OUTPUT_FASTA:
                        show_ARRAY(NUC)(stdout,calls,"",0);
                        break;
                }
                free_ARRAY(NUC)(rcseq);
            }
            unfiltered_count++;
        } else {
            // Represent filtered fasta's as N's
            if ( OUTPUT_FASTA==simopt->format ){
                const uint32_t lim = (model->paired)?(2*model->ncycle):model->ncycle;
                for ( uint32_t cy=0 ; cy<lim ; cy++){
                    show_NUC(stdout,NUC_AMBIG);
                }
            }
        }
        fputc('\n',stdout);
        if(NULL!=fpout){ fputc('\n',fpout); }
        free_SEQ(seq);
        seq_count++;
        if( (seq_count%1000)==0 ){ fprintf(stderr,"\rDone: %8u",seq_count); }
    }
    fprintf(stderr,"\rFinished generating %8u sequences\n",seq_count);
    if(simopt->purity_cycles>0){ fprintf(stderr,"%8u sequences passed filter.\n",unfiltered_count);}
    if(NULL!=fpout){fclose(fpout);}
    free_ARRAY(NUC)(calls);
    free_MAT(loglike);
    free_MAT(intensities);

    // Print error summary
    fputs("Summary of errors, calling by maximum likelihood\n",stderr);
    fputs("Cycle  Count  Phred   lower, upper",stderr);
    if(simopt->paired){ fputs("   Count  Phred   lower, upper",stderr);}
    for ( uint32_t i=0 ; i<model->ncycle ; i++){
        const real_t e = ((real_t)error[i])/unfiltered_count;
	    fprintf(stderr,"\n%3u: %7u %6.2f (%6.2f,%6.2f)",i+1,error[i], phred(e), phred(prop_upper(e,unfiltered_count)), phred(prop_lower(e,unfiltered_count)));
        if(simopt->paired){ 
            const real_t e2 = ((real_t)error2[i])/unfiltered_count;
            fprintf(stderr,"%7u %6.2f (%6.2f,%6.2f)",error2[i], phred(e2), phred(prop_upper(e2,unfiltered_count)), phred(prop_lower(e2,unfiltered_count)));
        }
    }
    fputc('\n',stderr);
    // Histograms
    fputs("Number of errors per read",stderr);
    for ( uint32_t i=0 ; i<6 ; i++){
        fprintf(stderr,"\n%2u: %7u %6.2f",i,errorhist[i],(100.0*errorhist[i])/unfiltered_count);
        if(simopt->paired){
            fprintf(stderr,"\t %7u %6.2f",errorhist[i],(100.0*errorhist[i])/unfiltered_count);
        }
    }
    fprintf(stderr,"\n>5: %7u %6.2f",errorhist[6],(100.0*errorhist[6])/unfiltered_count);
    if(simopt->paired){
        fprintf(stderr,"\t %7u %6.2f",errorhist[6],(100.0*errorhist[6])/unfiltered_count);
    }
    fputc('\n',stderr);
    free_MODEL(model);
    free_SIMOPT(simopt);

    return EXIT_SUCCESS;
}

