#ifndef _LAPACK_H
#define _LAPACK_H

#ifdef linux
#define F77_NAME(A) A ## _
#else
#define F77_NAME(A) A
#endif

const char LAPACK_UPPER[] = {'U'};
const char LAPACK_LOWER[] = {'L'};
const char LAPACK_UNITTRI[] = {'U'};
const char LAPACK_NONUNITTRI[] = {'N'};

// Float functions
void F77_NAME(spotrf)( const char * uplo, const int * n, float * A, const int * lda, int * info );
void F77_NAME(strtri)( const char * uplo, const char * diag, const int * n, float * a, const int * lda, int * info);

// Double functions
void F77_NAME(dpotrf)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dtrtri)( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);


// Generic definitions
#ifdef USEFLOAT
    #define potrf   F77_NAME(spotrf)
    #define trtri   F77_NAME(strtri)
#else
    #define potrf   F77_NAME(dpotrf)
    #define trtri   F77_NAME(dtrtri)
#endif

#endif

