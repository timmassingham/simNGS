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

#include <assert.h>
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
#include "kumaraswamy.h"
#include "lambda_distribution.h"


#define DEFAULT_DUST_PROB 0.0
#define STRINGIFY(A) #A
#define ILLUMINA_ADAPTER "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT"
#define PROGNAME "simNGS"
#define PROGVERSION "1.5"

enum paired_type { PAIRED_TYPE_SINGLE=0, PAIRED_TYPE_CYCLE, PAIRED_TYPE_PAIRED };
char * paired_type_str[] = {"single","cycle","paired"};
enum outformat { OUTPUT_LIKE=0, OUTPUT_FASTA, OUTPUT_FASTQ };
const char * output_format_str[] = { "likelihood", "fasta", "fastq" };


ARRAY(NUC) ambigseq = {NULL,0};
ARRAY(PHREDCHAR) ambigphred = {NULL,0};

typedef struct {
    MAT intensities;
    MAT loglike;
    ARRAY(NUC) calls;
    ARRAY(PHREDCHAR) quals;
    bool pass_filter;
} * CALLED;

typedef struct{
    CSTRING name;
    bool paired;
    real_t lambda1, lambda2;
    ARRAY(NUC) seq, rcseq;
    MAT int1,int2;
} * SEQSTR;

MAT reverse_complement_MAT(const MAT mat);
CALLED reverse_complement_CALLED( const CALLED called);
void free_CALLED(CALLED called);

void free_SEQSTR( SEQSTR seqstr){
    if(NULL==seqstr){ return;}
    free_MAT(seqstr->int1);
    free_MAT(seqstr->int2);
    free_CSTRING(seqstr->name);
    free_ARRAY(NUC)(seqstr->seq);
    free_ARRAY(NUC)(seqstr->rcseq);
    free(seqstr);
}

void show_SEQSTR(FILE * fp, const SEQSTR seqstr){
    if(NULL==fp){return;}
    if(NULL==seqstr){ return;}
    fprintf(fp,"Intensities for %s\n",seqstr->name);
    fprintf(fp,"Forwards, brightness = %f\n",seqstr->lambda1);
    show_ARRAY(NUC)(fp,seqstr->seq,"",40);
    fputc('\n',fp);
    show_MAT(fp,seqstr->int1,4,8);
    if(seqstr->paired){
        fprintf(fp,"Backwards, brightness = %f\n",seqstr->lambda2);
        show_ARRAY(NUC)(fp,seqstr->rcseq,"",40);
        fputc('\n',fp);
        show_MAT(fp,seqstr->int2,4,8);
    }
}
        
#define X(A) A ## SEQSTR
#include "circbuff.def"
#undef X


void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"" PROGNAME "\"\n"
"Simulate likelihoods for Illumina data from fasta format files\n"
"\n"
"Usage:\n"
"\t" PROGNAME " [-a adapter] [-b shape1:scale1:shape2:scale2] [-c correlation] [-d] [-D prob]\n"
"\t       [-f nimpure:ncycle:threshold] [-F factor] [-g prob] [-i filename] [-I]\n"
"\t       [-j range:a:b] [-l lane] [-m] [-M matrix file] [-N noise file] [-n ncycle]\n"
"\t       [-o output_format] [-p option] [-P phasing file] [-q quantile] [-r mu] [-R]\n"
"\t       [-s seed] [-t tile] [-v factor ] runfile [seq.fa ... ]\n"
"\t" PROGNAME " --help\n"
"\t" PROGNAME " --licence\n"
"\t" PROGNAME " --license\n"
"\t" PROGNAME " --version\n"
"" PROGNAME " reads from stdin and writes to stdout. Messages and progess\n"
"indicators are written to stderr.\n"
"\n"
"Example:\n"
"\tcat sequences.fa | " PROGNAME " runfile > sequences.fq\n"
"\t" PROGNAME " runfile seq1.fa seq2.fa > sequences.fq\n"
,fp);
}

void fprint_licence (FILE * fp){
    validate(NULL!=fp,);
    fputs(
"  " PROGNAME " software for simulating likelihoods for next-gen sequencing machines\n"
#include "copyright.inc"
    ,fp);
}

void fprint_version(FILE * fp){
    validate(NULL!=fp,);
    fputs(
"  " PROGNAME " software for simulating likelihoods for next-gen sequencing machines\n"
"Version " PROGVERSION " (compiled: " __DATE__ " using " __VERSION__ ")\n"
, fp);
}


void fprint_help( FILE * fp){
    validate(NULL!=fp,);
    fputs(
/*
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/
"\n"
"-a, --adapter sequence [default: " ILLUMINA_ADAPTER "]\n"
"\tSequence to pad reads with if they are shorter than the number of\n"
"cycles required, reflecting the adpater sequence used for sample preparaton.\n"
"A null adapter (i.e. pad with ambiguity characters) can be specified by -a \"\"\n" 
"\n" 
"-b, --brightness shape1:scale1[:shape2:scale2] [default: as runfile]\n"
"\tShape and scale of cluster brightness distribution for each end of pair.\n"
"Parameters for second end (shape2 and scale2) are optional; those for the first\n"
"end will be repeated if necessary.\n"
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
"-D, --dust probability [default: no dust]\n"
"\tProbability of dust occurring on a particular cycle, resulting in\n"
"extremely bright observations in the second \"C\" channel. Cross-talk, phasing\n"
"and noise matrices must be specified. A value of 1e-5 is typical.\n"
"\n"
"-f, --filter nimpure:ncycle:threshold [default: no filtering]\n"
"\tUse purity filtering on generated intensities, allowing a maximum of\n"
"nimpure cyles in the first ncycles with a purity greater than threshold.\n"
"\n"
"-F, --final factor [default: see below]\n"
"\tVariance adjustment for final cycle, derived from input covariance\n"
"matrix by default if the number of cycles required is fewer than that described\n"
"in the runfile or 1 if equal.\n"
"\n"
"-g, --generalised, --generalized probability [default: set from mutation rate]\n"
"\tProbability of a generalised error, a mistaken bsae not due to\n"
"base calling error. If not given, the probability of a generalised error\n"
"is set to the base mutation rate (see the --mutate option).\n"
"\n"
"-i, --intensities filename [default: none]\n"
"\tWrite the processed intensities generated to \"filename\".\n"
"\n"
"-I, --illumina\n"
"\tProduce Illumina scaled quality values where required. Ascii\n"
"representation of quality value is ascii(quality+64) rather than the more\n"
"usual ascii(quality+33).\n"
"\n"
"-j, --jumble range:a:b [default: none]\n"
"\tJumble generated intensities with those of other reads, simulating cases\n"
"where clusters have merged.\n"
"A random read is picked from the <range> reads after each each and their\n"
"intensities mixed to final intensities for analysis. Mixing is according to a\n"
"Kumaraswamy distribution with parameters a and b, a being the shape parameter\n"
"for the current read and b that for the randomly picked read.\n"
"\n"
"-l, --lane lane [default: as runfile]\n"
"\tSet lane number\n"
"\n"
"-m, --mutate, --mutate=insertion:deletion:mutation [default: 1e-5:1e-6:1e-4]\n"
"\tSimple model of sequence mutation to reflect sample preparation errors.\n"
"When the -m or --mutate options are given without an argument, the mutational\n"
"process is turned off otherwise the default parameters are used.\n"
"An alternative process of mutation may be specified using the format:\n"
"\t--mutate=1e-5:1e-6:1e-4\n"
"\n"
"-M, --matrix filename [default: none]\n"
"\tFile to read cross-talk matrix from. Not required for general\n"
"simulation of sequence and qualities.\n"
"\n"
"-n, --ncycles ncycles [default: as runfile]\n"
"\tNumber of cycles to do, up to maximum allowed for runfile.\n"
"\n"
"-N, --noise filename [default: none]\n"
"\tFile to read systematic noise matrix from. Not required for general\n"
"simulation of sequence and qualities.\n"
"\n"
"-o, --output format [default: likelihood]\n"
"\tFormat in which to output results. Either \"likelihood\", \"fasta\",\n"
"or \"fastq\".\n"
"\n"
"-p, --paired option [default: single]\n"
"\tTreat run as paired-end.\n"
"Valid options are: \"single\", \"paired\", \"cycle\".\n"
"For \"cycle\" reads, the results are reported in machine cycle order\n"
"(the second end is reverse complemented and appended to first).\n"
"The \"paired\" option splits the two ends of the read into separate\n"
"records, indicated by a suffix /1 or /2 to the read name. The second\n"
"end is reported in the opposite orientation to the first.\n"
"For single-ended runs treated as paired, the covariance matrix is\n"
"duplicated to make two uncorrelated pairs.\n"
"\n"
"-P, --phasing filename [default: none]\n"
"\tFile to read phasing matrix from. Not required for general\n"
"simulation of sequence and qualities.\n"
"\n"
"-q, --quantile quantile [default: 0]\n"
"\tQuantile below which cluster brightness is discarded and redrawn from\n"
"distribution.\n"
"\n"
"-r, --robust mu [default: 1e-6]\n"
"\tCalculate robustified likelihood, equivalent to adding mu to every\n"
"likelihood.\n"
"\n"
"-R, --raw\n"
"\tDump raw intensities rather than processed intensities.\n"
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
    { "dust",       required_argument, NULL, 'D' },
    { "final",      required_argument, NULL, 'F' },
    { "filter",     required_argument, NULL, 'f' },
    { "generalised", required_argument, NULL, 'g' },
    { "generalized", required_argument, NULL, 'g' },
    { "intensities", required_argument, NULL, 'i'},
    { "illumina",   no_argument,       NULL, 'I' },
    { "jumble",     required_argument, NULL, 'j' },
    { "lane",       required_argument, NULL, 'l' },
    { "mutate",     optional_argument, NULL, 'm' },
    { "matrix",	    required_argument, NULL, 'M' },
    { "ncycle",     required_argument, NULL, 'n' },
    { "noise",      required_argument, NULL, 'N' },
    { "output",     required_argument, NULL, 'o' },
    { "paired",     required_argument, NULL, 'p' },
    { "phasing",    required_argument, NULL, 'P' },
    { "quantile",   required_argument, NULL, 'q' },
    { "robust",     required_argument, NULL, 'r' },
    { "raw",        no_argument,       NULL, 'R' },
    { "seed",       required_argument, NULL, 's' },
    { "tile",       required_argument, NULL, 't' },
    { "variance",   required_argument, NULL, 'v' },
    { "help",       no_argument,       NULL, 'h' },
    { "licence",    no_argument,       NULL, 0 },
    { "license",    no_argument,       NULL, 0 },
    { "version",    no_argument,       NULL, 1 },
    { NULL, 0, NULL, 0}
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
    real_t corr,threshold;
    Distribution dist1,dist2;
    enum paired_type paired;
    bool desc;
    real_t mu,generr,sdfact;
    real_t final_factor[4];
    uint32_t seed;
    uint32_t tile,lane;
    real_t purity_threshold;
    uint32_t purity_cycles,purity_max;
    CSTRING intensity_fn;
    enum outformat format;
    bool mutate;
    real_t ins,del,mut;
    bool jumble;
    uint32_t bufflen;
    real_t a,b;
    ARRAY(NUC) adapter;
    real_t dustProb;
    MAT M,P,N;
    MAT invM,invP,Mt,Pt;
    bool illumina,dumpRaw;
} * SIMOPT;

SIMOPT new_SIMOPT(void){
    SIMOPT opt = calloc(1,sizeof(*opt));
    validate(NULL!=opt,NULL);
    // Set defaults
    opt->ncycle = 0;
    opt->dist1 = NULL;
    opt->dist2 = NULL;
    opt->corr = 1.0;
    opt->threshold = 0.0;
    opt->paired = PAIRED_TYPE_SINGLE;
    opt->desc = false;
    opt->mu = 1e-6;
    opt->generr= -1.0;
    opt->seed = 0;
    opt->sdfact = 1.0;
    opt->final_factor[0] = -1.0;
    opt->tile = 0;
    opt->lane = 0;
    opt->purity_threshold = 0.;
    opt->purity_cycles = 0;
    opt->purity_max = 0;
    opt->intensity_fn = NULL;
    opt->format = OUTPUT_FASTQ;
    opt->mutate = true;
    opt->ins=1e-5; opt->del=1e-6; opt->mut=1e-4;
    opt->jumble = false;
    opt->bufflen = 1; opt->a=0.; opt->b=0;
    opt->adapter = nucs_from_string(ILLUMINA_ADAPTER);
    opt->dustProb = 0.0;
    opt->M = opt->P = opt->N = NULL;
    opt->invM = opt->invP = NULL;
    opt->Mt = opt->Pt = NULL;
    opt->illumina = false;
    opt->dumpRaw = false;
    return opt;
}

void free_SIMOPT(SIMOPT opt){
    validate(NULL!=opt,);
    free(opt->intensity_fn);
    free_ARRAY(NUC)(opt->adapter);
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
    fprintf( fp,"paired\t%s\n",paired_type_str[simopt->paired]);
    fprintf( fp,"Brightness correlation\t%f\n",simopt->corr);
    fprintf( fp,"mu\t%f\n",simopt->mu);
    fputs("Lambda end1:\t",fp); show_Distribution(fp,simopt->dist1);
    fputs("Lambda end1:\t",fp); show_Distribution(fp,simopt->dist2);
    fprintf( fp,"threshold\t%f\n",simopt->threshold);
    fprintf( fp,"variance factor\t%f\n",simopt->sdfact*simopt->sdfact);
    fprintf( fp,"tile\t%u\tlane%u\n",simopt->tile,simopt->lane);
    fprintf( fp,"seed\t%u\n",simopt->seed);

    if(simopt->purity_cycles!=0){
       fprintf( fp,"Purity filtering: threshold %f. Maximum of %u inpure in %u cycles\n",simopt->purity_threshold, simopt->purity_max,simopt->purity_cycles);
    } else {
       fputs("No purity filtering.\n",fp);
    }

    if(simopt->mutate){
        fprintf(fp,"Input sequence will be mutated with ins %f, del %f, mut %f\n",simopt->ins,simopt->del,simopt->mut);
    }
    if(simopt->jumble){
        fprintf(fp,"Intensities will be jumbled with shape parameters %f , %f (range %u)\n",simopt->a,simopt->b,simopt->bufflen);
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
    
    while ((ch = getopt_long(argc, argv, "a:b:c:dD:F:f:g:i:Ij:l:mM:n:N:o:p:P:q:r:Rs:t:uv:h", longopts, NULL)) != -1){
        int ret;
        unsigned long int i=0,j=0;
	real_t param1[2],param2[2];
        switch(ch){
        case 'a':   free_ARRAY(NUC)(simopt->adapter);
                    simopt->adapter = nucs_from_string(optarg);
                    break;
	case 'b':   ret = sscanf(optarg,real_format_str ":" real_format_str ":" real_format_str ":" real_format_str,&param1[0],&param1[1],&param2[0],&param2[1]);
                    if(ret!=2 && ret!=4){ errx(EXIT_FAILURE,"Incorrect number of arguments for brightness (got %d).",ret);}
		    if(ret==2){
			    param2[0] = param1[0];
			    param2[1] = param1[1];
		    }
                    if(param1[0]<=0. || param2[0]<=0.){ errx(EXIT_FAILURE,"Brightness shape must be greater than zero."); }
                    if(param1[1]<=0. || param2[1]<=0.){ errx(EXIT_FAILURE,"Brightness scale must be greater than zero."); }
		    simopt->dist1 = new_Distribution('W',param1);
		    simopt->dist2 = new_Distribution('W',param2);
		    show_Distribution(stdout,simopt->dist1);
                    break;
        case 'c':   sscanf(optarg,real_format_str,&simopt->corr);
                    if(simopt->corr<-1.0 || simopt->corr>1.0){errx(EXIT_FAILURE,"Correlation between end brightness should be in [-1,1]. Was given %f.",simopt->corr);}
                    break;
        case 'd':   simopt->desc = true;
                    break;
	case 'D':   sscanf(optarg,real_format_str,&simopt->dustProb);
		    if(simopt->dustProb<0.0 || simopt->dustProb>1.0){errx(EXIT_FAILURE,"Dust probability should be in [0,1]. Was given %f.",simopt->dustProb);}
		    break;
	case 'F':   ret = sscanf(optarg, real_format_str ":" real_format_str ":" real_format_str ":" real_format_str ,&simopt->final_factor[0],&simopt->final_factor[1],&simopt->final_factor[2],&simopt->final_factor[3]);
		    if(ret!=4 && ret!=1){ errx(EXIT_FAILURE,"Incorrect number of arguments for final variance."); }
		    if(1==ret){
			    // If only read one double, copy over
			    simopt->final_factor[1] = simopt->final_factor[2] = simopt->final_factor[3] = simopt->final_factor[0];
		    }
		    if( simopt->final_factor[0]<=0.0 ||
		        simopt->final_factor[1]<=0.0 ||
			simopt->final_factor[2]<=0.0 ||
			simopt->final_factor[3]<=0.0 ){
		            errx(EXIT_FAILURE,"Arguments for final variance muct be positive.");
		    }
		    break;
        case 'f':   ret = sscanf(optarg, "%lu:%lu:" real_format_str,&i,&j,&simopt->purity_threshold);
                    if(ret!=3){ errx(EXIT_FAILURE,"Insufficient arguments for filtering.");}
                    simopt->purity_max = i;
                    simopt->purity_cycles = j;
                    if( simopt->purity_threshold<0. || simopt->purity_threshold>1.0){
                        errx(EXIT_FAILURE,"Purity threshold is %f but should be between 0 and 1.",simopt->purity_threshold);
                    }
                    break;
	case 'g':   simopt->generr = parse_real(optarg);
                    if(simopt->generr<0.0 || simopt->generr>1.0){errx(EXIT_FAILURE,"Generalized error is %f but must be a probability [0,1]",simopt->generr);}
		    break;
        case 'i':   simopt->intensity_fn = copy_CSTRING(optarg);
                    break;
                case 'j':   ret = sscanf(optarg, "%u:" real_format_str ":" real_format_str, &simopt->bufflen,&simopt->a, &simopt->b);
                    if( ret!=3 ){ errx(EXIT_FAILURE,"Insufficient arguments for jumbling.");}
                    if(0==simopt->bufflen){
                        errx(EXIT_FAILURE,"Range for jumbling must be positive");
                    }
                    if(0==simopt->a || 0==simopt->b){
                        errx(EXIT_FAILURE,"Jumbling not defined when shape parameters are zero.");
                    }
                    simopt->jumble = true;
                    break;
	case 'I':   simopt->illumina = true;
		    break;
        case 'l':   simopt->lane = parse_uint(optarg);
                    if(simopt->lane==0){errx(EXIT_FAILURE,"Lane number must be greater than zero.");}
                    break;
        case 'm':   if(NULL==optarg){ // No optional argument
			    simopt->mutate = false;
			    simopt->ins = simopt->del = simopt->mut = 0.0;
		    } // Optional argument present
		    else {
		    	ret = sscanf(optarg, real_format_str ":" real_format_str ":" real_format_str,&simopt->ins,&simopt->del,&simopt->mut);
                    	if( ret!=3 ){ errx(EXIT_FAILURE,"Insufficient arguments for mutation.");}
                    	if(!isprob(simopt->ins) || !isprob(simopt->del) || !isprob(simopt->mut) ){
                        	errx(EXIT_FAILURE,"Mutation parameters not probabilities. Given: ins %f, del %f, mut %f",simopt->ins,simopt->del,simopt->mut);
                    	}
                    	if(simopt->ins+simopt->del+simopt->mut>1.0){
                        	errx(EXIT_FAILURE,"Mutation parameters sum to greater than one.");
                    	}
                    	simopt->mutate = true;
		    }
                    break;
        case 'M':   simopt->M = new_MAT_from_file(optarg,0,0);
		    if(NULL==simopt->M){ errx(EXIT_FAILURE,"Failed to read cross-talk matrix from file %s",optarg); }
		    break;
        case 'n':   sscanf(optarg,"%u",&simopt->ncycle);
                    if(simopt->ncycle==0){errx(EXIT_FAILURE,"Number of cycles to simulate must be greater than zero.");}
                    break;
        case 'N':   simopt->N = new_MAT_from_file(optarg,0,0);
		    if(NULL==simopt->N){ errx(EXIT_FAILURE,"Failed to read noise matrix from file %s",optarg); }
		    break;
        case 'o':   if( strcasecmp(optarg,output_format_str[OUTPUT_LIKE])==0 ){ simopt->format = OUTPUT_LIKE; }
                    else if ( strcasecmp(optarg,output_format_str[OUTPUT_FASTA])==0 ){ simopt->format = OUTPUT_FASTA; }
                    else if ( strcasecmp(optarg,output_format_str[OUTPUT_FASTQ])==0 ){
                        simopt->format = OUTPUT_FASTQ;
                        if(simopt->mu==0){ simopt->mu = 1e-5;}
                    } else {
                        errx(EXIT_FAILURE,"Unrecognised output option %s.",optarg);
                    }
                    break;
        case 'p':   if( strcasecmp(optarg,paired_type_str[PAIRED_TYPE_SINGLE])==0 ){
                        simopt->paired = PAIRED_TYPE_SINGLE;
                    } else if ( strcasecmp(optarg,paired_type_str[PAIRED_TYPE_PAIRED])==0 ){
                        simopt->paired = PAIRED_TYPE_PAIRED;
                    } else if ( strcasecmp(optarg,paired_type_str[PAIRED_TYPE_CYCLE])==0 ){
                        simopt->paired = PAIRED_TYPE_CYCLE;
                    } else {
                        errx(EXIT_FAILURE,"Unrecognised paired option %s.",optarg);
                    }
                    break;
        case 'P':   simopt->P = new_MAT_from_file(optarg,0,0);
		    if(NULL==simopt->P){ errx(EXIT_FAILURE,"Failed to read phasing matrix from file %s",optarg); }
	            break;
        case 'q':   simopt->threshold = parse_real(optarg);
                    if(!isprob(simopt->threshold) ){ 
                       errx(EXIT_FAILURE,"Threshold quantile to discard brightness must be a probability (got %e)\n",simopt->threshold);
                    }
                    break;
        case 'r':   simopt->mu = parse_real(optarg);
                    if(simopt->mu<0.0){errx(EXIT_FAILURE,"Robustness \"mu\" must be non-negative.");}
                    break;
	case 'R':   simopt->dumpRaw = true;
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
        case 1:
            fprint_version(stderr);
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

void output_likelihood_sub(const SIMOPT simopt, const uint32_t x, const uint32_t y, const CALLED called1, const CALLED called2){
    const bool has_called2 = (called2!=NULL);
    fprintf(stdout,"%u\t%u\t%u\t%u",simopt->lane,simopt->tile,x,y);
    if(called1->pass_filter){
        fprint_intensities(stdout,"",called1->loglike,false);
        if(has_called2){fprint_intensities(stdout,"",called2->loglike,false);}
    }
    fputc('\n',stdout);
}

void output_likelihood(const SIMOPT simopt, const uint32_t x, const uint32_t y, const CALLED called1, const CALLED called2){
	switch(simopt->paired){
	case PAIRED_TYPE_SINGLE:
	case PAIRED_TYPE_CYCLE:
		output_likelihood_sub(simopt,x,y,called1,called2);
		break;
	case PAIRED_TYPE_PAIRED:
		output_likelihood_sub(simopt,x,y,called1,NULL);
		output_likelihood_sub(simopt,x,y,called2,NULL);
		break;
	default:
		errx(EXIT_FAILURE,"Unrecognised case %s (%s:%d)",__func__,__FILE__,__LINE__);
	}	
}

void output_fasta_sub(const SIMOPT simopt, const char * seqname, const char * suffix, const CALLED called1, const CALLED called2){
    const bool has_called2 = (called2!=NULL);
    fprintf(stdout,">%s%s\n",seqname,suffix);
    if(called1->pass_filter){
        show_ARRAY(NUC)(stdout,called1->calls,"",0);
        if(has_called2){ show_ARRAY(NUC)(stdout,called2->calls,"",0);}
    } else {
        show_ARRAY(NUC)(stdout,ambigseq,"",0);
        if(has_called2){ show_ARRAY(NUC)(stdout,ambigseq,"",0);}
    }
    fputc('\n',stdout);
}

void output_fastq_sub(const SIMOPT simopt, const char * seqname, const char * suffix, const CALLED called1, const CALLED called2){
    const bool has_called2 = (called2!=NULL);
    fprintf(stdout,"@%s%s\n",seqname,suffix);
    if(called1->pass_filter){
        show_ARRAY(NUC)(stdout,called1->calls,"",0);
        if(has_called2){ show_ARRAY(NUC)(stdout,called2->calls,"",0);}
    } else {
        show_ARRAY(NUC)(stdout,ambigseq,"",0);
        if(has_called2){ show_ARRAY(NUC)(stdout,ambigseq,"",0);}
    }
    fputs("\n+\n",stdout);
    if(called1->pass_filter){
        show_ARRAY(PHREDCHAR)(stdout,called1->quals,"",0);
        if(has_called2){ show_ARRAY(PHREDCHAR)(stdout,called2->quals,"",0);}
    } else {
        show_ARRAY(PHREDCHAR)(stdout,ambigphred,"",0);
        if(has_called2){ show_ARRAY(PHREDCHAR)(stdout,ambigphred,"",0);}
    }
    fputc('\n',stdout);
}

void output_fasta(const SIMOPT simopt, const char * seqname, const CALLED called1, const CALLED called2){
        switch(simopt->paired){
        case PAIRED_TYPE_SINGLE:
        case PAIRED_TYPE_CYCLE:
                output_fasta_sub(simopt,seqname,"",called1,called2);
                break;
        case PAIRED_TYPE_PAIRED:
                output_fasta_sub(simopt,seqname,"/1",called1,NULL);
                output_fasta_sub(simopt,seqname,"/2",called2,NULL);
                break;
        default:
                errx(EXIT_FAILURE,"Unrecognised case %s (%s:%d)",__func__,__FILE__,__LINE__);
        }
}



void output_fastq(const SIMOPT simopt, const char * seqname, const CALLED called1, const CALLED called2){
	switch(simopt->paired){
	case PAIRED_TYPE_SINGLE:
        case PAIRED_TYPE_CYCLE:
		output_fastq_sub(simopt,seqname,"",called1,called2);
		break;
	case PAIRED_TYPE_PAIRED:
		output_fastq_sub(simopt,seqname,"/1",called1,NULL);
		output_fastq_sub(simopt,seqname,"/2",called2,NULL);
		break;
	default:
		errx(EXIT_FAILURE,"Unrecognised case %s (%s:%d)",__func__,__FILE__,__LINE__);
	}
}

void output_results(FILE * intout, const SIMOPT simopt, const char * seqname, const uint32_t x, const uint32_t y, const CALLED called1, const CALLED called2){
    // Output raw intensities if required
    if(NULL!=intout){
        fprintf(intout,"%u\t%u\t%u\t%u",simopt->lane,simopt->tile,x,y);
        fprint_intensities(intout,"",called1->intensities,false);
        if(NULL!=called2){fprint_intensities(intout,"",called2->intensities,false);}
        fputc('\n',intout);
    }

    // Output in format requested
    switch(simopt->format){
       case OUTPUT_LIKE:
           output_likelihood(simopt,x,y,called1,called2);
           break;
       case OUTPUT_FASTA:
           output_fasta(simopt,seqname,called1,called2);
           break;
       case OUTPUT_FASTQ:
           output_fastq(simopt,seqname,called1,called2);
           break;
       default:
           errx(EXIT_FAILURE,"Unrecognised format in %s (%s:%d)",__func__,__FILE__,__LINE__);
    }
}

struct pair_double { double x1,x2;};

struct pair_double correlated_distribution(const real_t threshold, const real_t corr, const Distribution dist1, const Distribution dist2){

    // Pick lambda using Gaussian Copula
    real_t lambda1=NAN,lambda2=NAN;
    real_t px=0.0,py=0.0;
    do{
        // Two correlated Gaussians
        real_t x = rstdnorm();
        real_t y = corr*x + sqrt(1-corr*corr) * rstdnorm();
        // Convert to uniform deviates (the copula)
        px = pstdnorm(x,false,false);
        py = pstdnorm(y,false,false);
    } while(px<threshold || py<threshold);
    // Convert to observation via inversion formula
    lambda1 = qdistribution(px,dist1,false,false);
    lambda2 = qdistribution(py,dist2,false,false);
    return (struct pair_double){lambda1,lambda2};
}


MAT mix_intensities(const MAT int1, const MAT int2, const real_t prop){
    if(NULL==int1 || NULL==int2){ return NULL;}
    validate(int1->nrow==int2->nrow && int1->ncol==int2->ncol,NULL);
    validate(isprob(prop),NULL);

    MAT intmix = new_MAT(int1->nrow,int1->ncol);
    if(NULL==intmix){ return NULL;}
    
    const uint32_t nelt = int1->nrow * int1->ncol;
    for ( uint32_t i=0 ; i<nelt ; i++){
        intmix->x[i] = int2->x[i] + prop * (int1->x[i]-int2->x[i]);
    }
    return intmix;
}

MAT random_mix_intensities(const MAT int1, const MAT int2, const real_t shape1, const real_t shape2){
    if(NULL==int1 || NULL==int2){ return NULL;}
    real_t prop = rkumaraswamy(shape1,shape2);
    return mix_intensities(int1,int2,prop);
}

    
CALLED process_intensities( MAT intensities, const real_t lambda, const MAT * invchol, const SIMOPT simopt){
    CALLED cl = calloc(1,sizeof(*cl));
    if(NULL==cl){ return NULL;}
    cl->intensities = intensities;
    cl->loglike = likelihood_cycle_intensities(simopt->sdfact,simopt->mu,lambda,intensities,invchol,NULL);
    cl->calls = call_by_maximum_likelihood(cl->loglike,cl->calls);
    cl->quals = quality_from_likelihood(cl->loglike,cl->calls,simopt->generr,simopt->illumina,cl->quals);
    cl->pass_filter = number_inpure_cycles(intensities,simopt->purity_threshold,simopt->purity_cycles) <= simopt->purity_max;
    if(simopt->dumpRaw){
        cl->intensities = unprocess_intensities(cl->intensities,simopt->Mt, simopt->Pt, simopt->N, NULL);
	free_MAT(intensities);
    }
    return cl;
}

MAT reverse_complement_MAT(const MAT mat){
	if(NULL==mat){ return NULL;}
	if(mat->nrow!=NBASE){
		errx(EXIT_FAILURE,"reverse_complement_MAT applied to incompatable matrix of %" SCNu32 "x%" SCNu32, mat->nrow,mat->ncol);
	}
	MAT m = new_MAT(mat->nrow,mat->ncol);
	if(NULL==m){ return NULL;}
	for ( uint32_t col=0 ; col<mat->ncol ; col++){
		const uint32_t idx = col*NBASE;
		const uint32_t idx_new = (mat->ncol - col - 1)*NBASE;
		m->x[idx_new+NUC_A] = mat->x[idx+NUC_T];
		m->x[idx_new+NUC_C] = mat->x[idx+NUC_G];
		m->x[idx_new+NUC_G] = mat->x[idx+NUC_C];
		m->x[idx_new+NUC_T] = mat->x[idx+NUC_A];
	}
	return m;
}

CALLED reverse_complement_CALLED( const CALLED called){
	if(NULL==called){ return NULL;}
	ARRAY(NUC) rcnuc = null_ARRAY(NUC);
	ARRAY(PHREDCHAR) revqual = null_ARRAY(PHREDCHAR);
	MAT rcintensities=NULL,rcloglike=NULL;
	CALLED called_new = NULL;

	// Convert entries of previous CALLED object
	rcnuc = reverse_complement(called->calls);
	if(rcnuc.elt==NULL){ goto cleanup;}
	revqual = reverse_quality(called->quals);
	if(revqual.elt==NULL){ goto cleanup;}
	rcintensities = reverse_complement_MAT(called->intensities);
	if(rcintensities==NULL){ goto cleanup;}
	rcloglike = reverse_complement_MAT(called->loglike);
	if(rcloglike==NULL){ goto cleanup;}
	
	// Create new CALLED object
	called_new = calloc(1,sizeof(*called_new));
	if(NULL==called_new){ goto cleanup;}
	called_new->intensities = rcintensities;
	called_new->loglike = rcloglike;
	called_new->calls = rcnuc;
	called_new->quals = revqual;
	called_new->pass_filter = called->pass_filter;

	return called_new;

cleanup:
	free(called_new);
	free_MAT(rcloglike);
	free_MAT(rcintensities);
	free_ARRAY(PHREDCHAR)(revqual);
	free_ARRAY(NUC)(rcnuc);
	return NULL;
}

void free_CALLED(CALLED called){
    if(NULL==called){ return;}
    free_MAT(called->loglike);
    free_MAT(called->intensities);
    free_ARRAY(NUC)(called->calls);
    free_ARRAY(PHREDCHAR)(called->quals);
    free(called);
}

void update_error_counts(const ARRAY(NUC) calls, const ARRAY(NUC) seq, uint32_t * error, uint32_t * errorhist){
    uint32_t ncycle = calls.nelt;
    uint32_t nerr = 0;
    for ( uint32_t i=0 ; i<ncycle ; i++){
        if(i>=seq.nelt || calls.elt[i] != seq.elt[i]){ nerr++; error[i]++;}
    }
    errorhist[(nerr<6)?nerr:6]++;
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
    argc--;
    argv++;
    if( simopt->desc ){
        show_MODEL(stderr,model);
        return EXIT_SUCCESS;
    } 
    fprintf(stderr,"Description of runfile:\n%s",model->label);
    
    // Resolve options and model
    // Should factor out into separate routine
    if(NULL!=simopt->dist1){
	    free_Distribution(model->dist1);
	    model->dist1 = copy_Distribution(simopt->dist1);
	    free_Distribution(model->dist2);
	    model->dist2 = copy_Distribution(simopt->dist2);
    }
    if(simopt->generr==-1.0){
	    simopt->generr = (simopt->mutate)?simopt->mut:0.0;
    }

    // Dust simulation requires phasing, cross-talk and noise matrices
    if(0.0!=simopt->dustProb){
        if(NULL==simopt->M || NULL==simopt->P || NULL==simopt->N){
	    errx(EXIT_FAILURE,"Cross-talk, phasing and noise matrices required to simulate dust");
	}
    }
    if(simopt->dumpRaw){
        if(NULL==simopt->M || NULL==simopt->P || NULL==simopt->N){
            errx(EXIT_FAILURE,"Cross-talk, phasing and noise matrices required to dump raw intensities");
	}
    }
    // Check that M, P and N dimensions are consistent with run file
    if(NULL!=simopt->M){
        if(NBASE!=simopt->M->nrow || NBASE!=simopt->M->ncol){
            errx(EXIT_FAILURE,"Cross-talk matrix has wrong dimension, got %d,%d",simopt->M->nrow,simopt->M->ncol);
	}
	// Want inverse for calculations
	simopt->invM = invert_MAT(simopt->M);
	simopt->Mt = transpose(simopt->M);
    }
    if(NULL!=simopt->P){
        if(model->ncycle!=simopt->P->nrow || model->ncycle!=simopt->P->ncol){
	    errx(EXIT_FAILURE,"Phasing matrix has wrong dimension, got %d,%d",simopt->P->nrow,simopt->P->ncol);
	}
	// Want inverse for calculations
	simopt->invP = invert_MAT(simopt->P);
	simopt->Pt = transpose(simopt->P);
    }
    if(NULL!=simopt->N){
	if(NBASE!=simopt->N->nrow || model->ncycle!=simopt->N->ncol){
	    errx(EXIT_FAILURE,"Systematic noise matrix has wrong dimension, got %d,%d",simopt->N->nrow,simopt->N->ncol);
	}
    }


    
    if(model->paired && simopt->paired==PAIRED_TYPE_SINGLE){
        fputs("Treating paired-end model as single-ended.\n",stderr);
        model->paired = false;
        free_MAT(model->cov2);
        model->cov2 = NULL;
    } else if(!model->paired && simopt->paired!=PAIRED_TYPE_SINGLE){
        fputs("Treating single-ended model as paired-end.\n",stderr);
        model->paired = true;
        model->cov2 = copy_MAT(model->cov1);
        model->chol2 = copy_MAT(model->chol1);
        model->invchol2 = calloc(model->ncycle,sizeof(*model->invchol2));
	if(!model->dist2){
		model->dist2 = copy_Distribution(model->dist1);
	}
        for ( uint32_t i=0 ; i<model->ncycle ; i++){
            model->invchol2[i] = copy_MAT(model->invchol1[i]);
        }
    }

    if(simopt->ncycle==0){
	    simopt->ncycle = model->ncycle;
    }

    if(simopt->ncycle>model->ncycle){
        fprintf(stderr,"Asked for more cycles than runfile allows. Doing %u.\n",model->ncycle);
    } else {
        MODEL newmodel = trim_MODEL(simopt->ncycle,simopt->final_factor,model);
        free_MODEL(model);
        model = newmodel; 
    }

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
    SEQ seq = NULL;
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
    
    // Create sequence of ambiguities for filtered calls
    ambigseq = new_ARRAY(NUC)(model->ncycle);
    ambigphred = new_ARRAY(PHREDCHAR)(model->ncycle);
    for( uint32_t i=0 ; i<model->ncycle ; i++){
        ambigseq.elt[i] = NUC_AMBIG;
        ambigphred.elt[i] = '!';
    }

    // Circular buffer for intensities. Size one if no buffer.
    CIRCBUFF(SEQSTR) circbuff = new_circbuff_SEQSTR(simopt->bufflen);
    FILE * fp = stdin;
    do { // Iterate through filenames
        if(argc>0){
            fp = fopen(argv[0],"r");
            if(NULL==fp){
                warnx("Failed to open file \"%s\" for input",argv[0]);
            }
        }
        while (NULL!=fp && (seq=sequence_from_fasta(fp))!=NULL){
            //show_SEQ(stderr,seq);
            if(simopt->mutate){
                SEQ mut = mutate_SEQ(seq,simopt->ins,simopt->del,simopt->mut);
                free_SEQ(seq);
                seq = mut;
            }

            SEQSTR seqstr = calloc(1,sizeof(*seqstr));
            seqstr->name = copy_CSTRING(seq->name);
            seqstr->seq = copy_ARRAY(NUC)(seq->seq);
            seqstr->paired = model->paired;
            free_SEQ(seq); seq=NULL;
            // Pick copula
            struct pair_double lambda = correlated_distribution(simopt->threshold,simopt->corr,model->dist1,model->dist2);
            seqstr->lambda1 = lambda.x1;
            seqstr->lambda2 = lambda.x2;

            // Generate intensities
            seqstr->int1 = generate_pure_intensities(simopt->sdfact,lambda.x1,seqstr->seq,simopt->adapter,model->ncycle,model->chol1_cycle,simopt->dustProb,simopt->invM,simopt->invP,simopt->N,NULL);
            if ( model->paired ){
                seqstr->rcseq = reverse_complement(seqstr->seq);
                seqstr->int2 = generate_pure_intensities(simopt->sdfact,lambda.x2,seqstr->rcseq,simopt->adapter,model->ncycle,model->chol2_cycle,simopt->dustProb,simopt->invM,simopt->invP,simopt->N,NULL);
            }
            // Store in buffer
            SEQSTR popped = push_circbuff_SEQSTR(circbuff,seqstr);
            if( NULL!=popped ){
                MAT intensities=NULL,intensities2=NULL;
                CALLED called1=NULL, called2=NULL;

                // Can only pop when buffer is full
                if( simopt->jumble ){
                    real_t prop = rkumaraswamy(simopt->a,simopt->b);
                    uint32_t randelt = (uint32_t)(circbuff->maxelt*runif());
                    assert(randelt>=0 && randelt<circbuff->maxelt);
                    intensities  = mix_intensities(popped->int1,circbuff->elt[randelt]->int1,prop);
                    intensities2 = mix_intensities(popped->int2,circbuff->elt[randelt]->int2,prop);
                } else {
                    intensities  = copy_MAT(popped->int1);
                    intensities2 = copy_MAT(popped->int2);
                }
            
                called1 = process_intensities(intensities,popped->lambda1,model->invchol1,simopt);
                update_error_counts(called1->calls,popped->seq,error,errorhist);
            
                called2 = process_intensities(intensities2,popped->lambda2,model->invchol2,simopt);
                update_error_counts(called2->calls,popped->rcseq,error2,errorhist2);
            
                if(called1->pass_filter){ unfiltered_count++;}
                uint32_t x = (uint32_t)( 1794 * runif());
                uint32_t y = (uint32_t)( 2048 * runif());
                output_results(fpout,simopt,popped->name,x,y,called1,called2);
            
                free_CALLED(called1); called1=NULL; intensities=NULL;
                free_CALLED(called2); called2=NULL; intensities2=NULL;
                free_SEQSTR(popped);

                seq_count++;
                if( (seq_count%1000)==0 ){ fprintf(stderr,"\rDone: %8u",seq_count); }        
            }
        }
        fclose(fp);
        argc--;
        argv++;
    } while(argc>0);
    // Buffer still contains (upto) simopt->bufflen elements Output.
    {
        MAT intensities=NULL,intensities2=NULL;
        CALLED called1=NULL, called2=NULL;
        uint32_t maxelt = (circbuff->maxelt<circbuff->nseen)?circbuff->maxelt:circbuff->nseen;
        uint32_t oldest = (circbuff->maxelt<circbuff->nseen)?(circbuff->nseen%circbuff->maxelt):0;
        for ( uint32_t i=0 ; i<maxelt ; i++){
            uint32_t idx = (i+oldest)%circbuff->maxelt;
            SEQSTR popped = circbuff->elt[idx];
            if( simopt->jumble ){
                real_t prop = rkumaraswamy(simopt->a,simopt->b);
                uint32_t randelt = idx;
                do {
                    randelt = (uint32_t)(maxelt*runif());
                } while(randelt==idx);
                assert(randelt>=0 && randelt<maxelt);
                intensities  = mix_intensities(popped->int1,circbuff->elt[randelt]->int1,prop);
                intensities2 = mix_intensities(popped->int2,circbuff->elt[randelt]->int2,prop);
            } else {
                intensities  = copy_MAT(popped->int1);
                intensities2 = copy_MAT(popped->int2);
            }

            called1 = process_intensities(intensities,popped->lambda1,model->invchol1,simopt);
            update_error_counts(called1->calls,popped->seq,error,errorhist);

            called2 = process_intensities(intensities2,popped->lambda2,model->invchol2,simopt);
            update_error_counts(called2->calls,popped->rcseq,error2,errorhist2);

            if(called1->pass_filter){ unfiltered_count++;}
            uint32_t x = (uint32_t)( 1794 * runif());
            uint32_t y = (uint32_t)( 2048 * runif());
            output_results(fpout,simopt,popped->name,x,y,called1,called2);

            free_CALLED(called1); called1=NULL; intensities=NULL;
            free_CALLED(called2); called2=NULL; intensities2=NULL;

            seq_count++;
            if( (seq_count%1000)==0 ){ fprintf(stderr,"\rDone: %8u",seq_count); }
        }
    }
    // Empty and free buffer
    for ( uint32_t i=0 ; i<circbuff->maxelt ; i++){
        free_SEQSTR(circbuff->elt[i]);
    }
    free_circbuff_SEQSTR(circbuff);
    
    fprintf(stderr,"\rFinished generating %8u sequences\n",seq_count);
    if(simopt->purity_cycles>0){ fprintf(stderr,"%8u sequences passed filter.\n",unfiltered_count);}
    if(NULL!=fpout){fclose(fpout);}
    free_ARRAY(PHREDCHAR)(ambigphred);
    free_ARRAY(NUC)(ambigseq);

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
    free(error); free(error2);
    // Histograms
    fputs("Number of errors per read",stderr);
    for ( uint32_t i=0 ; i<6 ; i++){
        fprintf(stderr,"\n%2u: %7u %6.2f%%",i,errorhist[i],(100.0*errorhist[i])/unfiltered_count);
        if(simopt->paired){
            fprintf(stderr,"\t %7u %6.2f%%",errorhist2[i],(100.0*errorhist2[i])/unfiltered_count);
        }
    }
    fprintf(stderr,"\n>5: %7u %6.2f",errorhist[6],(100.0*errorhist[6])/unfiltered_count);
    if(simopt->paired){
        fprintf(stderr,"\t %7u %6.2f",errorhist2[6],(100.0*errorhist2[6])/unfiltered_count);
    }
    fputc('\n',stderr);
    free_MODEL(model);
    free_SIMOPT(simopt);

    return EXIT_SUCCESS;
}

