/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the simNGS software for simulating likelihoods
 *  for next-generation sequencing machines.
 *
 *  AYB is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AYB is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with simNGS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _LAPACK_H
#define _LAPACK_H

#ifdef __linux
#define F77_NAME(A) A ## _
#else
#define F77_NAME(A) A
#endif

static const char LAPACK_UPPER[] = {'U'};
static const char LAPACK_LOWER[] = {'L'};
static const char LAPACK_UNITTRI[] = {'U'};
static const char LAPACK_NONUNITTRI[] = {'N'};
static const char LAPACK_TRANS[] = {'T'};
static const char LAPACK_NOTRANS[] = {'N'};
static const int LAPACK_UNIT[] = {1};

// Float functions
void F77_NAME(spotrf)( const char * uplo, const int * n, float * A, const int * lda, int * info );
void F77_NAME(strtri)( const char * uplo, const char * diag, const int * n, float * a, const int * lda, int * info);
void F77_NAME(strmv)( const char * uplo, const char * trans, const char * diag, const int * n, const float * A, const int * lda, float * x, const int * incx);
void F77_NAME(sposv)(const char * uplo, const int * N, const int * NRHS, float * A, const int * lda, float * B, const int * ldb, const int * info);
void F77_NAME(sposvx)(const char * fact, const char * uplo, const int * N,
                      const int * NRHS, float * A, const int * lda, 
                      const int * af, const int * ldaf, char * equed,
                      float * S, float * B, const int * ldb, float * x,
                      const int * ldx, float * rcond, float * ferr,
                      float * berr, float * work, int * iwork, int * info);
void F77_NAME(sgelss)( const int * M, const int * N, const int * NRHS, float * A,
                       const int * lda, float * B, const int * ldb, float * S,
                       float * RCOND, int * rank, float * WORK, int * LWORK, int * INFO);
void F77_NAME(sgemm)(   const char * transa, const char * transb, 
                        const int * M, const int * N, const int *K,
                        const float * alpha, const float * A, const int * lda,
                        const float * B, const int * ldb, const float * beta,
                        float * C, const int * ldc);
void F77_NAME(sgemv)(   const char * trans, const int * M, const int * N,
                        const float * alpha, const float * A, const int * lda,
                        const float * X, const int * incx, const float * beta,
                        float * y, const int * incy);
void F77_NAME(sgetrf)(const int * M, const int * N, float * A, const int * lda,
                 int * ipiv, int * info);

void F77_NAME(sgetri)(const int * N, float * A, const int * lda, const int * ipiv,
                 float * work, const int * lwork, int * info);


// Double functions
void F77_NAME(dpotrf)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dtrtri)( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);
void F77_NAME(dtrmv)( const char * uplo, const char * trans, const char * diag, const int * n, const double * A, const int * lda, double * x, const int * incx);
void F77_NAME(dposv)(const char * uplo, const int * N, const int * NRHS, double * A, const int * lda, double * B, const int * ldb, const int * info);
void F77_NAME(dposvx)(const char * fact, const char * uplo, const int * N,
                      const int * NRHS, double * A, const int * lda, 
                      const int * af, const int * ldaf, char * equed,
                      double * S, double * B, const int * ldb, double * x,
                      const int * ldx, double * rcond, double * ferr,
                      double * berr, double * work, int * iwork, int * info);
void F77_NAME(dgelss)( const int * M, const int * N, const int * NRHS, double * A,
                       const int * lda, double * B, const int * ldb, double * S,
                       double * RCOND, int * rank, double * WORK, int * LWORK, int * INFO);
void F77_NAME(dgemm)(   const char * transa, const char * transb, 
                        const int * M, const int * N, const int *K,
                        const double * alpha, const double * A, const int * lda,
                        const double * B, const int * ldb, const double * beta,
                        double * C, const int * ldc);
void F77_NAME(dgemv)(   const char * trans, const int * M, const int * N,
                        const double * alpha, const double * A, const int * lda,
                        const double * X, const int * incx, const double * beta,
                        double * y, const int * incy);
void F77_NAME(dgetrf)(const int * M, const int * N, double * A, const int * lda,
                 int * ipiv, int * info);

void F77_NAME(dgetri)(const int * N, double * A, const int * lda, const int * ipiv,
                 double * work, const int * lwork, int * info);



// Generic definitions
#ifdef USEFLOAT
    #define potrf   F77_NAME(spotrf)
    #define trtri   F77_NAME(strtri)
    #define trmv    F77_NAME(strmv)
    #define posv    F77_NAME(sposv)
    #define posvx   F77_NAME(sposvx)
    #define gelss   F77_NAME(sgelss)
    #define gemm    F77_NAME(sgemm)
    #define gemv    F77_NAME(sgemv)
    #define getrf   F77_NAME(sgetrf)
    #define getri   F77_NAME(sgetri)
#else
    #define potrf   F77_NAME(dpotrf)
    #define trtri   F77_NAME(dtrtri)
    #define trmv    F77_NAME(dtrmv)
    #define posv    F77_NAME(dposv)
    #define posvx   F77_NAME(dposvx)
    #define gelss   F77_NAME(dgelss)
    #define gemm    F77_NAME(dgemm)
    #define gemv    F77_NAME(dgemv)
    #define getrf   F77_NAME(dgetrf)
    #define getri   F77_NAME(dgetri)
#endif

#endif

