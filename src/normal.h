#ifndef _NORMAL_H
#define _NORMAL_H

#ifndef _UTILITY_H
#include "utility.h"
#endif

#ifndef _STDINT_H
#include <stdint.h>
#endif

#ifndef _MATRIX_H
#include "matrix.h"
#endif

real_t rstdnorm( void );
real_t dstdnorm( const real_t x, const bool logd);

real_t rnormal( const real_t mean, const real_t sd );
real_t dnorm( const real_t x, const real_t m, const real_t sd, const bool logd);

MAT rmultinorm( const MAT mean, const MAT At, const uint32_t n, MAT z);
real_t dmultinorm( const MAT x, const MAT mean, const MAT At, const uint32_t n, const bool logd);

#endif

