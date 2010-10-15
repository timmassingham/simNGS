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
#include <tgmath.h>
#include <ctype.h>
#include "sequence.h"
#include "random.h"
#include "normal.h"

enum strand_opt { STRAND_RANDOM, STRAND_SAME, STRAND_OPPOSITE };

#define Q_(A) #A
#define QUOTE(A) Q_(A)
#define DEFAULT_COV         0.055
#define DEFAULT_INSERT      400
#define DEFAULT_NCYCLE      45
#define DEFAULT_COVERAGE    2.0
#define DEFAULT_BIAS        0.5
#define PROGNAME "simLibrary"
#define PROGVERSION "1.2"

uint32_t nfragment_from_coverage(const uint32_t genlen, const real_t coverage, const uint32_t readlen, const bool paired){
    const uint32_t bases_per_read = paired?(2*readlen):readlen;
    return (uint32_t)(0.5+(genlen*coverage)/bases_per_read);
}

CSTRING fragname(const CSTRING name,const char strand, const uint32_t loc, const uint32_t fraglen){
    char * str;
    asprintf(&str,"%s (Strand %c Offset %u--%u)",name,strand,loc,loc+fraglen);
    return str;
}


void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"" PROGNAME "\"\n"
"Split sequence into a simulated library of fragments\n"
"\n"
"Usage:\n"
"\t" PROGNAME " [-b bias] [-c cov] [-i insertlen] [-m multiplier_file]\n"
"\t         [-n nfragments] -p [-r readlen] [-s strand] [-v variance]\n"
"\t         [-x coverage] [--seed seed] seq1.fa ...\n"
"\t" PROGNAME " --help\n"
"\t" PROGNAME " --licence\n"
"\t" PROGNAME " --version\n"
PROGNAME " reads from stdin and writes to stdout. Messages and progess\n"
"indicators are written to stderr.\n"
"\n"
"Example:\n"
"\tcat genome.fa | " PROGNAME " > library.fa\n"
,fp);
}

void fprint_licence(FILE * fp){
    validate(NULL!=fp,);
    fputs(
"  " PROGNAME " software for simulating libraries of fragments from genomic sequence\n"
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
"-b, --bias bias [default: " QUOTE(DEFAULT_BIAS) "]\n"
"\tStrand bias for sampling. The probability of sampling a read from the\n"
"positive strand.\n"
"\n"
"-c, --cov cov [default: " QUOTE(DEFAULT_COV) "]\n"
"\tCoefficient Of Variance (COV) for the read lengths. The COV is the\n"
"ratio of the variance to the mean^2 and is related to the variance of the\n"
"log-normal distribution by cov = exp(var)-1. If the variance option is set,\n"
"it takes presidence.\n"
"\n"
"-i, --insert insert_length [default: " QUOTE(DEFAULT_INSERT) "]\n"
"\tMean length of insert. The mean length of the reads sampled is the\n"
"insert length plus twice the read length.\n"
"\n"
"-m, --multipliers file [default: see text]\n"
"\tA file containing multiplers, one per line, for the number of\n"
"fragments to produce for each input sequence. The total number of\n"
"fragments produced for each sequence is controlled by the --coverage\n"
"and --nfragments but scaled by the multiplier. If there are multiple\n"
"input sequences and a multipliers file is not given, or contains\n"
"insufficient multipliers, then a default multiplier of one is used.\n"
"\n"
"-n, --nfragments nfragments [default: from coverage]\n"
"\tNumber of fragments to produce for library. By default the number of\n"
"fragments is sufficient for the coverage given. If the number of fragments\n"
"is set then this option takes priority.\n"
"\n"
"-p, --paired [default: true ]\n"
"\tTurn off paired-end generation. The average fragment length will be\n"
"shorter by an amount equal to the read length but the main effect of turning\n"
"off paired-end generation is in the coverage calculations: twice as many\n"
"fragments will be generated for single-ended runs as paired-end.\n"
"\n"
"-r, --readlen read_length [default: " QUOTE(DEFAULT_NCYCLE) "]\n"
"\tRead length to sample. Affects the total length of fragments produced\n"
"and the total number of fragments produced via the coverage.\n"
"\n"
"-s, --strand strand [default: random]\n"
"\tStrand from which simulated reads are to come from, relative\n"
"to current strand. Options are: opposite, random, same.\n"
"\n"
"--seed seed [default: clock]\n"
"\tSet seed from random number generator.\n"
"\n"
"-v, --variance variance [default: from COV]\n"
"\tThe variance of the read length produced. By default, the variance is\n"
"set using the effective read length and the Coefficient of Variance so the\n"
"standard deviation is proportional to the mean. Setting the variance takes\n"
"priority over the COV.\n"
"\n"
"-x, --coverage coverage [ default: " QUOTE(DEFAULT_COVERAGE) "]\n"
"\tAverage coverage of original sequence for simulated by fragments. If\n"
"the number of fragments to produce is set, it take priority over the coverage.\n"
,fp);
}

static struct option longopts[] = {
    { "bias",       required_argument, NULL, 'b'},
    { "cov",        required_argument, NULL, 'c'},
    { "insert",     required_argument, NULL, 'i'},
    { "multipliers",required_argument, NULL, 'm'},
    { "nfragments", required_argument, NULL, 'n'},
    { "paired",     no_argument,       NULL, 'p'},
    { "readlen",    required_argument, NULL, 'r'},
    { "strand",     required_argument, NULL, 's'},
    { "seed",       required_argument, NULL, 2 },
    { "variance",   required_argument, NULL, 'v'},
    { "coverage",   required_argument, NULL, 'x'},
    { "help",       no_argument,       NULL, 'h'},
    { "licence",    no_argument,       NULL, 0 },
    { "version",    no_argument,       NULL, 1 }
};

typedef struct {
    uint32_t insertlen, ncycle, nfragment;
    bool paired;
    uint32_t seed;
    real_t variance,cov,strand_bias,coverage;
    enum strand_opt strand;
    FILE * multiplier_fp;
} * OPT;

OPT new_OPT(void){
    OPT opt = calloc(1,sizeof(*opt));
    validate(NULL!=opt,NULL);
    opt->insertlen = DEFAULT_INSERT;
    opt->ncycle = DEFAULT_NCYCLE;
    opt->seed = 0;
    opt->variance = 0;
    opt->cov = DEFAULT_COV;
    opt->coverage = DEFAULT_COVERAGE;
    opt->nfragment = 0;
    opt->paired = true;
    opt->strand_bias = DEFAULT_BIAS;
    opt->strand = STRAND_RANDOM;
    opt->multiplier_fp = NULL;
    return opt;
}

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

OPT parse_options(const int argc, char * const argv[] ){
    int ch;
    OPT opt = new_OPT();
    validate(NULL!=opt,NULL);
    
    while ((ch = getopt_long(argc, argv, "b:c:i:m:n:pr:s:v:x:h", longopts, NULL)) != -1){
        switch(ch){
        case 'b':
            opt->strand_bias = parse_real(optarg);
            if(!isprob(opt->strand_bias)){ errx(EXIT_FAILURE,"Positive strand bias should be between zero and one, got %f",opt->strand_bias);}
            break;
        case 'c':
            opt->cov = parse_real(optarg);
            if(opt->cov<0.0){errx(EXIT_FAILURE,"Coefficient Of Variance of for insert size should be non-zero");}
            break;
        case 'i':
            opt->insertlen = parse_uint(optarg);
            if(0==opt->insertlen){errx(EXIT_FAILURE,"Insert length should be strictly positive");}
            break;
        case 'm':
            opt->multiplier_fp = fopen(optarg,"r");
            if(NULL==opt->multiplier_fp){
                warnx("Failed to open multiplier file \"%s\". Using multiplier of one.");
            }
            break;
        case 'n':
            opt->nfragment = parse_uint(optarg);
            break;
        case 'p':
            opt->paired = false;
            break;
        case 'r':
            opt->ncycle = parse_uint(optarg);
            break;
        case 's':
            if(!strncasecmp(optarg,"same",4)){ opt->strand = STRAND_SAME; break;}
            if(!strncasecmp(optarg,"opposite",8)){ opt->strand = STRAND_OPPOSITE; break;}
            if(!strncasecmp(optarg,"random",6)){ opt->strand = STRAND_RANDOM; break;}
            errx(EXIT_FAILURE,"Unrecognised choice of strand \"%s\"",optarg); break;
        case 2:
            opt->seed = parse_uint(optarg);
            break;
        case 'v':
            opt->variance = parse_real(optarg);
            if(opt->variance<0.0){errx(EXIT_FAILURE,"Variance of insert size should be non-zero");}
            break;
        case 'x':
            opt->coverage = parse_real(optarg);
            if(0==opt->coverage){errx(EXIT_FAILURE,"Coverage should be strictly positive");}
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
    return opt;
}

int main ( int argc, char * argv[]){
    
    OPT opt = parse_options(argc,argv);
    if(NULL==opt){
        errx(EXIT_FAILURE,"Failed to parse options");
    }
    argc -= optind;
    argv += optind;
    
    if ( opt->seed==0 ){
        uint32_t seed = (uint32_t) time(NULL);
        fprintf(stderr,"Using seed %u\n",seed);
        opt->seed = seed;
    }
    init_gen_rand( opt->seed );

    // Alter strand_bias if complete bias is required
    switch(opt->strand){
    case STRAND_SAME: opt->strand_bias = 1.0; break;
    case STRAND_OPPOSITE: opt->strand_bias = 0.0; break;
    default: break;
    }
    
    real_t effectivelen = opt->insertlen + opt->ncycle + ((opt->paired)?opt->ncycle:0);
    real_t log_sd = (0==opt->variance)?
        sqrt(log1p(opt->cov)) :
        sqrt(log1p(opt->variance/(effectivelen*effectivelen)));
    real_t log_mean = log(effectivelen) - 0.5 * log_sd * log_sd;
    
    FILE * fp = stdin;
    SEQ seq = NULL;
    uint32_t skipped_seq = 0, tot_fragments = 0;
    do { // Iterate through filenames
        if(argc>0){
            fp = fopen(argv[0],"r");
            if(NULL==fp){
                warnx("Failed to open file \"%s\" for input",argv[0]);
            }
        }
        while (NULL!=fp && (seq=sequence_from_fasta(fp))!=NULL){
            SEQ rcseq = reverse_complement_SEQ(seq);
            // Read multiplier from file, if available
            double multiplier = 1.;
            if(NULL!=opt->multiplier_fp){
                int ret = fscanf(opt->multiplier_fp,"%lf",&multiplier);
                // Clean-up if file has run out of multipliers.
                if(1!=ret){
                    warnx("Failed to read multiplier from file. Will use default of 1 from now on"); 
                    fclose(opt->multiplier_fp);
                    opt->multiplier_fp = NULL;
                }
                if(multiplier<0.0){
                    warnx("Invalid fragment multiplier %lf. Using 1");
                    multiplier = 1.0;
                }
            }

            // Number of fragments
            uint32_t nfragment = multiplier * ( (opt->nfragment)?opt->nfragment:nfragment_from_coverage(seq->length,opt->coverage,opt->ncycle,opt->paired) );
            for ( uint32_t i=0 ; i<nfragment ; i++,tot_fragments++){
                const uint32_t fraglen = (opt->paired)?(uint32_t)(exp(rnorm(log_mean,log_sd))):opt->ncycle;
                if(fraglen>seq->length){
                     if(0==skipped_seq){
                         warnx("Length of fragment (2*readlen+insert) is greater than sequence length. Skipping");
                     }
                     skipped_seq++;
                     continue;
                }
                const uint32_t loc = (uint32_t)((seq->length-fraglen)*runif()); // Location is uniform
                char strand = (runif()<opt->strand_bias)?'+':'-';

                SEQ sampseq = (strand=='+')?sub_SEQ(seq,loc,fraglen):sub_SEQ(rcseq,loc,fraglen);
                CSTRING sampname = fragname(seq->name,strand,loc,fraglen);
                free_CSTRING(sampseq->name);
                sampseq->name = sampname;
                show_SEQ(stdout,sampseq);
                free_SEQ(sampseq);
                if( (tot_fragments%100000)==99999 ){ fprintf(stderr,"\rDone: %8u",tot_fragments+1); }
            }
            free_SEQ(rcseq);
            free_SEQ(seq);
        }
        fclose(fp);
        argc--;
        argv++;
    } while(argc>0);
    fprintf(stderr,"\rFinished %8u\n",tot_fragments);
    if(skipped_seq>0){
        fprintf(stderr,"Skipped %" SCNu32 " fragments.\n",skipped_seq);
    }

    
    return EXIT_SUCCESS;
}
