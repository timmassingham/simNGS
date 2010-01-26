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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "utility.h"

/* Basic functions for basic generic type */
#define _FREE(_TYPE)    void free_ ## _TYPE  ( _TYPE _val){ return; }
#define _COPY(_TYPE)    _TYPE copy_ ## _TYPE ( const _TYPE _val ){ return _val; }

/* char */
_FREE(char)
_COPY(char)
char new_char( void ){ return 0;}
void show_char( FILE * fp, const char c){ fputc(c,fp); }

/* bool */
_FREE(bool)
_COPY(bool)
bool new_bool( void ){ return 0;}
void show_bool( FILE * fp, const bool c){ fprintf(fp,(c==true)?"true":"false"); }

/* int */
_FREE(int)
_COPY(int)
int new_int( void ){ return 0;}
void show_int( FILE * fp, const int i){ fprintf(fp,"%d",i); }

/* int8_t */
_FREE(int8_t)
_COPY(int8_t)
int8_t new_int8_t( void ){ return 0;}
void show_int8_t( FILE * fp, const int8_t i){ fprintf(fp,"%d",i); }

/* int16_t */
_FREE(int16_t)
_COPY(int16_t)
int16_t new_int16_t( void ){ return 0;}
void show_int16_t( FILE * fp, const int16_t i){ fprintf(fp,"%d",i); }

/* int32_t */
_FREE(int32_t)
_COPY(int32_t)
int32_t new_int32_t( void ){ return 0;}
void show_int32_t( FILE * fp, const int32_t i){ fprintf(fp,"%d",i); }

/* int64_t */
_FREE(int64_t)
_COPY(int64_t)
int64_t new_int64_t( void ){ return 0;}
void show_int64_t( FILE * fp, const int64_t i){ fprintf(fp,"%d",i); }

/* uint8_t */
_FREE(uint8_t)
_COPY(uint8_t)
uint8_t new_uint8_t( void ){ return 0;}
void show_uint8_t( FILE * fp, const uint8_t i){ fprintf(fp,"%u",i); }

/* uint16_t */
_FREE(uint16_t)
_COPY(uint16_t)
uint16_t new_uint16_t( void ){ return 0;}
void show_uint16_t( FILE * fp, const uint16_t i){ fprintf(fp,"%u",i); }

/* uint32_t */
_FREE(uint32_t)
_COPY(uint32_t)
uint32_t new_uint32_t( void ){ return 0;}
void show_uint32_t( FILE * fp, const uint32_t i){ fprintf(fp,"%u",i); }

/* uint64_t */
_FREE(uint64_t)
_COPY(uint64_t)
uint64_t new_uint64_t( void ){ return 0;}
void show_uint64_t( FILE * fp, const uint64_t i){ fprintf(fp,"%u",i); }

/* float */
_FREE(float)
_COPY(float)
float new_float(void){ return 0.;}
void show_float( FILE * fp, const float f){ fprintf(fp,"%f",f);}

/* double */
_FREE(double)
_COPY(double)
double new_double(void){ return 0.;}
void show_double( FILE * fp, const double f){ fprintf(fp,"%f",f);}

/* real_t */
_FREE(real_t)
_COPY(real_t)
real_t new_real_t(void){ return 0.;}
void show_real_t( FILE * fp, const real_t f){ fprintf(fp,"%f",f);}

/* cstring */
CSTRING new_CSTRING(const size_t len){ return calloc(len+1,sizeof(char)); }
void free_CSTRING( CSTRING c){ safe_free(c); }
CSTRING copy_CSTRING( const CSTRING c){
    validate(NULL!=c,NULL);
    CSTRING n = new_CSTRING( strlen(c) );
    validate(NULL!=n,NULL);
    strcpy(n,c);
    return n;
}
void show_CSTRING(FILE * fp, const CSTRING c){
    validate(NULL!=fp,);
    validate(NULL!=c,);
    fputs(c,fp);
}
