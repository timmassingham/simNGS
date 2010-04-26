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
 

#include <tgmath.h>
#include <stdbool.h>
#include "random.h"
#include "utility.h"

static real_t pow1p( real_t x, real_t y){
    return exp(y*log1p(x));
}

real_t pkumaraswamy(  real_t x,  real_t shape1,  real_t shape2,  bool tail,  bool logp){
    validate(shape1>0,NAN);
    validate(shape2>0,NAN);

    real_t s = -expm1(shape1*log(x));
    real_t log_tail_p = shape2 * log(s);
    // Tail calculations
    if(true==tail){
        if(true==logp){ return log_tail_p;}
        else          { return exp(log_tail_p);}
    }
    // Normal calculations
    real_t prob = -expm1(log_tail_p);
    if(false==logp){
        return prob;
    }
    return log(prob);
}

real_t qkumaraswamy(  real_t p,  real_t shape1,  real_t shape2,  bool tail,  bool logp){
    validate(shape1>0,NAN);
    validate(shape2>0,NAN);
    
    real_t lp = NAN;
    if(false==tail && false==logp){ lp = pow1p(-p,1.0/shape1); }
    if(false==tail && true==logp) { lp = exp(log(-expm1(p))/shape1);}
    if(true==tail  && false==logp){ lp = pow(p,1.0/shape1);}
    if(true==tail  && true==logp) { lp = exp(p/shape1);}
        
    return pow1p(-lp,1.0/shape2);
}

real_t dkumaraswamy(  real_t x,  real_t shape1,  real_t shape2,  bool logd ){
    validate(shape1>0,NAN);
    validate(shape2>0,NAN);

    real_t s = -expm1(shape1*log(x));
    real_t log_dens = log(shape1) + log(shape2) + (shape1-1.0)*log(x) + (shape2-1.0)*log(s);
    return (logd) ? log_dens : exp(log_dens);
}

real_t rkumaraswamy(  real_t shape1,  real_t shape2 ){
    validate(shape1>0,NAN);
    validate(shape2>0,NAN);
    /* Uses property that U uniform implies 1-U uniform
     * on the assumption that log(p) will be faster that log1p
     * on most platforms (more effort spent optimising). See
     * qkumaraswamy.
     */
    return qkumaraswamy(runif(),shape1,shape2,true,false);
}

#ifdef TEST
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

bool hassuffix( const char * str, const char * suf ){
    const uint32_t suf_len = strlen(suf);
    const uint32_t str_len = strlen(str);
    if ( suf_len>str_len){ return false;}
    return (0==strcmp(str+(str_len-suf_len),suf));
}
    
int main ( int argc, char * argv[] ){
    if(1==argc){
        fputs("Usage:\n"
"qkumaraswamy x shape1 shape2 tail logp\n"
"pkumaraswamy p shape1 shape2 tail logp\n"
"dkumaraswamy x shape1 shape2 log\n"
"rkumaraswamy shape1 shape2 seed\n", stderr);
        return EXIT_FAILURE;
    }
    
    double x,p,shape1,shape2;
    unsigned int tail,logp;
    long unsigned int seed;
    printf("%s",argv[0]);

    if ( hassuffix(argv[0],"qkumaraswamy") ){
        sscanf(argv[1],"%le",&x);
        sscanf(argv[2],"%le",&shape1);
        sscanf(argv[3],"%le",&shape2);
        sscanf(argv[4],"%u",&tail);
        sscanf(argv[5],"%u",&logp);
        printf("(%e,%e,%e,%u,%u)=",x,shape1,shape2,tail,logp);
        printf("%e\n",qkumaraswamy(x,shape1,shape2,tail,logp));
        return EXIT_SUCCESS;
    }

    if ( hassuffix(argv[0],"pkumaraswamy") ){
        sscanf(argv[1],"%le",&p);
        sscanf(argv[2],"%le",&shape1);
        sscanf(argv[3],"%le",&shape2);
        sscanf(argv[4],"%u",&tail);
        sscanf(argv[5],"%u",&logp);
        printf("(%e,%e,%e,%u,%u)=",p,shape1,shape2,tail,logp);
        printf("%e\n",pkumaraswamy(p,shape1,shape2,tail,logp));
        return EXIT_SUCCESS;
    }
    
    if ( hassuffix(argv[0],"dkumaraswamy") ){
        sscanf(argv[1],"%le",&x);
        sscanf(argv[2],"%le",&shape1);
        sscanf(argv[3],"%le",&shape2);
        sscanf(argv[4],"%u",&logp);
        printf("(%e,%e,%e,%u)=",x,shape1,shape2,logp);
        printf("%e\n",dkumaraswamy(x,shape1,shape2,logp));
        return EXIT_SUCCESS;
    }
    
    if ( hassuffix(argv[0],"rkumaraswamy") ){
        sscanf(argv[1],"%le",&shape1);
        sscanf(argv[2],"%le",&shape2);
        sscanf(argv[3],"%lu",&seed);
        init_gen_rand(seed);
        printf("(%e,%e,%lu)=",shape1,shape2,seed);
        printf("%e\n",rkumaraswamy(shape1,shape2));
        return EXIT_SUCCESS;
    }
    
    return EXIT_FAILURE;
}
#endif

