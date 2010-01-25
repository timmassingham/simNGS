#ifndef _UTILITY_H
#define _UTILITY_H

#ifndef _STDLIB_H
#include <stdlib.h>
#endif
#ifndef _STDIO_H
#include <stdio.h>
#endif
#ifndef _STDINT_H
#include <stdint.h>
#endif
#ifndef _STDBOOL_H
#include <stdbool.h>
#endif
#ifndef _ERR_H
#include <err.h>
#endif
#ifndef _SYS_ERRNO_H_
#include <errno.h>
#endif


#define safe_free(A) { if(NULL==(A)){free(A);} (A)=NULL; }


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


#ifdef USEFLOAT
    typedef float real_t;
    #define real_format_str "%f"
#else
    typedef double real_t;
    #define real_format_str "%lf"
#endif

#define BASIC_INTERFACE(_TYPE)                  \
    _TYPE new_ ## _TYPE   ( void);              \
    void free_ ## _TYPE  ( _TYPE _val);         \
    _TYPE copy_ ## _TYPE ( const _TYPE _val);   \
    void show_ ## _TYPE  ( FILE * fp, const _TYPE _val)

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
void free_CSTRING(void);
CSTRING copy_CSTRING(const CSTRING cstr);
void show_CSTRING(FILE *fp, const CSTRING cstr);

#endif
