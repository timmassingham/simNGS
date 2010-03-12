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

// Double functions
void F77_NAME(dpotrf)( const char * uplo, const int * n, double * A, const int * lda, int * info );
void F77_NAME(dtrtri)( const char * uplo, const char * diag, const int * n, double * a, const int * lda, int * info);
void F77_NAME(dtrmv)( const char * uplo, const char * trans, const char * diag, const int * n, const double * A, const int * lda, double * x, const int * incx);


// Generic definitions
#ifdef USEFLOAT
    #define potrf   F77_NAME(spotrf)
    #define trtri   F77_NAME(strtri)
    #define trmv    F77_NAME(strmv)
#else
    #define potrf   F77_NAME(dpotrf)
    #define trtri   F77_NAME(dtrtri)
    #define trmv    F77_NAME(dtrmv)
#endif

#endif

