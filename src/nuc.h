#ifndef _NUC_H
#define _NUC_H

#include <stdio.h>
#include <stdint.h>

typedef char NUC;
#define NBASE       4

#define NUC_AMBIG   4
#define NUC_A       0
#define NUC_C       1
#define NUC_G       2
#define NUC_T       3

typedef char PHREDCHAR;

void show_NUC(FILE * fp, const NUC nuc);
    
NUC nuc_from_char( const char c);
char char_from_nuc(const NUC nuc);
NUC complement(const NUC nuc);
NUC * reverse_complement(const NUC * nuc, const uint32_t len);
PHREDCHAR phredchar_from_char( const char c);

#endif

