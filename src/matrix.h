/* Standard copyright header here */
/* Description of module */

#ifndef MATRIX_H
#define MATRIX_H

#ifndef _STDINT_H
#include <stdint.h>
#endif

#ifndef _UTILITY_H
#include "utility.h"
#endif

/*  Specification did not specify what type of input should be taken. Use a
 * typedef so it can be easily changed, although care must be taken with the
 * input and output routines (printf etc) to make sure types match.
 */ 

struct _matrix_str {
    int    nrow, ncol;
    real_t    * x;
    };

// Make future abstraction easier
typedef struct _matrix_str * MAT;

// create/free
MAT new_MAT( const int nrow, const int ncol );
void free_MAT( MAT mat );
MAT copy_MAT( const MAT mat);
MAT copyinto_MAT( MAT matout, const MAT matin);
MAT new_MAT_from_array( const uint32_t nrow, const uint32_t ncol, const real_t * x);

// Input, output
MAT new_MAT_from_file( const char * filename );
MAT new_MAT_from_fp(FILE * fp );
void fprint_MAT( FILE * fp, const MAT mat);
void show_MAT( FILE * fp, const MAT mat, const uint32_t mrow, const uint32_t mcol);

// Identities
bool is_square(const MAT mat);

// Special matrices
MAT identity_MAT( const int nrow);

// Operations
MAT vectranspose( const MAT mat, const unsigned int p );
MAT reshape_MAT( MAT mat, const int nrow);
MAT cholesky( MAT mat);
MAT invert_cholesky( MAT mat);
MAT trim_MAT( MAT mat, const int mrow, const int mcol, const bool forwards);
MAT * block_diagonal_MAT( const MAT mat, const int n);
MAT scale_MAT(MAT mat, const real_t f);
#endif
