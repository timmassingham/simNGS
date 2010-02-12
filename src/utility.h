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

#ifndef _UTILITY_H
#define _UTILITY_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <err.h>
#include <errno.h>


#ifdef USEFLOAT
    typedef float real_t;
    #define real_format_str "%f"
#else
    typedef double real_t;
    #define real_format_str "%lf"
#endif

bool isprob( const real_t p);

#define safe_free(A) { if(NULL!=(A)){free(A);} (A)=NULL; }
#define SWAP(A,B) { typeof(A) tmp = A; A=B; B=tmp; }


#ifdef FAILEARLY
    #define validate(A,B) { if( !(A) ){ fprintf(stderr,"Validation failure for %s in %s at %s:%d",#A,__func__,__FILE__,__LINE__); abort();} }
#elif defined(SKIPVALIDATE)
    #define validate(A,B)
#else
    #define validate(A,B) { if( !(A) ){ return B; } }
#endif

#ifndef _TIME_H
#include <time.h>
#endif

static void timestamp(const char * str, FILE * fp){
    time_t t = time(&t);
    char * c = ctime(&t);
    c[24]='\t';
    fprintf(fp,"%s%s",c,str);
}

static inline bool bool_from_int( const int n){
    if(n!=0 && n!=1){
        errx(EINVAL,"Invalid boolean %d",n);
    }
    if(n==0){ return false; }
    return true;
}




#define BASIC_INTERFACE(_TYPE)                  \
    _TYPE new_ ## _TYPE   ( void);              \
    void free_ ## _TYPE  ( _TYPE _val);         \
    _TYPE copy_ ## _TYPE ( const _TYPE _val);   \
    void show_ ## _TYPE  ( FILE * fp, const _TYPE _val); \
    _TYPE read_ ## _TYPE ( FILE *  fp);

BASIC_INTERFACE(char);
BASIC_INTERFACE(int);
BASIC_INTERFACE(bool);
// Signed ints
BASIC_INTERFACE(int8_t);
BASIC_INTERFACE(int16_t);
BASIC_INTERFACE(int32_t);
BASIC_INTERFACE(int64_t);
// Unsigned ints
BASIC_INTERFACE(uint8_t);
BASIC_INTERFACE(uint16_t);
BASIC_INTERFACE(uint32_t);
BASIC_INTERFACE(uint64_t);
// Float, doubles, reals
BASIC_INTERFACE(float);
BASIC_INTERFACE(double);
BASIC_INTERFACE(real_t);
// Simple strings
typedef char * CSTRING;
CSTRING new_CSTRING(const size_t len);
void free_CSTRING(CSTRING cstr);
CSTRING copy_CSTRING(const CSTRING cstr);
void show_CSTRING(FILE *fp, const CSTRING cstr);
CSTRING read_CSTRING(FILE *fp);

int skipUntilChar ( FILE * fp, const char c);

#endif
