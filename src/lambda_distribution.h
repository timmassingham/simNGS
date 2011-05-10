#ifndef LAMBDA_DISTRIBUTIONS_H
#define LAMBDA_DISTRIBUTIONS_H

#include "utility.h"
#include <stdbool.h>

int nparameter_distribution(const char dist);
bool validate_parameters(const real_t * param , const char dist);
real_t qdistribution(const real_t px, const char dist, const real_t * param, const bool tail, const bool logp);
real_t qlogistic(const real_t p, const real_t loc, const real_t sc, const bool tail, const bool logp);

#endif /* LAMBDA_DISTRIBUTIONS_H */
