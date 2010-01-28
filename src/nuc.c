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
#include "utility.h"
#include "nuc.h"

void show_NUC(FILE * fp, const NUC nuc){
    validate(NULL!=fp,);
    fputc(char_from_nuc(nuc),fp);
}
        
NUC nuc_from_char( const char c){
    switch(c){
        case 'A':
        case 'a': return NUC_A;
        case 'C':
        case 'c': return NUC_C;
        case 'G':
        case 'g': return NUC_G;
        case 'T':
        case 't': return NUC_T;
        case 'N':
        case 'n': return NUC_AMBIG;
        default:
            fprintf(stderr,"Unrecognised nucleotide \"%c\". Returning NUC_AMBIG\n",c);
    }
    return NUC_AMBIG;
}

char char_from_nuc(const NUC nuc){
    switch(nuc){
    case NUC_A: return 'A';
    case NUC_C: return 'C';
    case NUC_G: return 'G';
    case NUC_T: return 'T';
    }
    return 'N';
}
    

NUC complement(const NUC nuc){
    switch(nuc){
    case NUC_A: return NUC_T;
    case NUC_C: return NUC_G;
    case NUC_G: return NUC_C;
    case NUC_T: return NUC_A;
    }
    return NUC_AMBIG;
}

ARRAY(NUC) reverse_complement(const ARRAY(NUC) nucs){
    validate(NULL!=nucs.elt,NULL_ARRAY(NUC));
    ARRAY(NUC) new_nuc = NEW_ARRAY(NUC)(nucs.nelt);
    validate(NULL!=new_nuc.elt,new_nuc);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        new_nuc.elt[i] = complement(nucs.elt[nucs.nelt-i-1]);
    }
    return new_nuc;
}

PHREDCHAR phredchar_from_char( const char c){
    validate(c>32,33);
    validate(c<127,126);
    return c;
}


uint32_t nucs_in_flows( const ARRAY(FLOW) flows){
    uint32_t n=0;
    for ( uint32_t i=0 ; i<flows.nelt ; i++){
        n += flows.elt[i];
    }
    return n;
}

uint32_t flows_in_nucs( const ARRAY(NUC) nucs, const ARRAY(NUC) flow_order ){
    validate(nucs.nelt==flow_order.nelt,0);
    uint32_t flow=0;
    for(uint32_t i=0 ; i<nucs.nelt && flow<flow_order.nelt; i++){
        while(flow<flow_order.nelt && nucs.elt[i]!=flow_order.elt[flow]){ flow++;}
    }
    return flow;
}


ARRAY(NUC) nucs_from_flow( const ARRAY(FLOW) flows, const ARRAY(NUC) flow_order){
    validate(flows.nelt==flow_order.nelt,NULL_ARRAY(NUC));
    uint32_t len = nucs_in_flows(flows);
    ARRAY(NUC) nucs = NEW_ARRAY(NUC)(len);
    if(NULL==nucs.elt){ return nucs;}
    for ( uint32_t flow=0,idx=0 ; flow<flows.nelt && idx<len ; flow++ ){
        for ( uint32_t j=0 ; j<flows.elt[flow] && idx<len ; j++,idx++){
            nucs.elt[idx] = flow_order.elt[flow];
        }
    }
    return nucs;
}

ARRAY(FLOW) flows_from_nucs( const ARRAY(NUC) nucs, const ARRAY(NUC) flow_order){
    ARRAY(FLOW) flows = NEW_ARRAY(FLOW)(flow_order.nelt);
    uint32_t flow=0;
    
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        while(flow<flow_order.nelt && nucs.elt[i]!=flow_order.elt[flow]){ flow++;}
        if(flow<flow_order.nelt){ flows.elt[flow]++;}
    }
    // Remaining elements of flows are zero by definition
    return flows;
}

#ifdef TEST
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main ( int argc, char * argv[]){
    if(argc!=3){
        fputs("Usage: test sequence flows\n",stderr);
        return EXIT_FAILURE;
    }
    uint32_t nnuc = strlen(argv[1]);
    ARRAY(NUC) nucs = NEW_ARRAY(NUC)(nnuc);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        nucs.elt[i] = nuc_from_char(argv[1][i]);
    }
    uint32_t nflow = strlen(argv[2]);
    ARRAY(NUC) flow_order = NEW_ARRAY(NUC)(nflow);
    for ( uint32_t i=0 ; i<nflow ; i++){
        flow_order.elt[i] = nuc_from_char(argv[2][i]);
    }
    
    fputs("Read sequence:     ",stdout);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        fputc(char_from_nuc(nucs.elt[i]),stdout);
    }
    fputc('\n',stdout);
    
    ARRAY(NUC) rc = reverse_complement(nucs);
    fputs("Reversed sequence: ",stdout);
    for ( uint32_t i=0 ; i<rc.nelt ; i++){
        fputc(char_from_nuc(rc.elt[i]),stdout);
    }
    fputc('\n',stdout);
    
    ARRAY(FLOW) flow = flows_from_nucs(nucs,flow_order);
    fputs("Sequence as flows: ",stdout);
    for ( uint32_t i=0 ; i<flow.nelt ; i++){
        fprintf(stdout,"%d",flow.elt[i]);
    }
    fputc('\n',stdout);
    
    ARRAY(NUC) trans = nucs_from_flow(flow,flow_order);
    fputs("Back-translated  : ",stdout);
    for ( uint32_t i=0 ; i<trans.nelt ; i++){
        fputc(char_from_nuc(trans.elt[i]),stdout);
    }
    fputc('\n',stdout);
}
#endif

