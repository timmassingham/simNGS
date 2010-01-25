#ifndef _LAPACK_H
#define _LAPACK_H

const char LAPACK_UPPER[] = {'U'};
const char LAPACK_LOWER[] = {'L'};
const char LAPACK_UNITTRI[] = {'U'};
const char LAPACK_NONUNITTRI[] = {'N'};

// Float functions
void spotrf( const char * uplo, const int * n, float * A, const int * lda, int * info );
void strtri( const char * uplo, const char * diag, const int * n, float * a, const int * lda, int * info);

// Double functions
void dpotrf( const char * uplo, const int * n, double * A, const int * lda, int * info );
void dtrtri( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);


// Generic definitions
#ifdef USEFLOAT
    #define potrf   spotrf
    #define trtri   strtri
#else
    #define potrf   dpotrf
    #define trtri   dtrtri
#endif

#endif

