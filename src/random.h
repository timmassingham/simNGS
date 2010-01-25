#ifndef _RANDOM_H
#define _RANDOM_H

#ifndef SFMT_H
#include "SFMT-src-1.3/SFMT.h"
#endif

#ifndef _UTILITY_H
#include "utility.h"
#endif

//Uniform RV on (0,1)                       
inline static real_t runif(void){
	return (((double)gen_rand64()) + 0.5) * (1.0/18446744073709551616.0L);
}

#endif

