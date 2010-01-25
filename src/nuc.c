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

NUC * reverse_complement(const NUC * nuc, const uint32_t len){
    validate(NULL!=nuc,NULL);
    NUC * new_nuc = calloc(len,sizeof(NUC));
    validate(NULL!=new_nuc,NULL);
    for ( uint32_t i=0 ; i<len ; i++){
        new_nuc[i] = complement(nuc[len-i-1]);
    }
    return new_nuc;
}

PHREDCHAR phredchar_from_char( const char c){
    validate(c>32,33);
    validate(c<127,126);
    return c;
}