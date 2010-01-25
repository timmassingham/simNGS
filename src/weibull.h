#ifndef _WEIBULL_H
#define _WEIBULL_H

#ifndef _STDBOOL_H
#include <stdbool.h>
#endif

#ifndef _UTILITY_H
#include "utility.h"
#endif

real_t pweibull( const real_t x, const real_t shape, const real_t scale, const bool tail, const bool logp);
real_t qweibull( const real_t p, const real_t shape, const real_t scale, const bool tail, const bool logp);
real_t dweibull( const real_t x, const real_t shape, const real_t scale, const bool logd );
real_t rweibull( const real_t shape, const real_t scale );

#endif
