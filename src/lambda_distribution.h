#ifndef LAMBDA_DISTRIBUTIONS_H
#define LAMBDA_DISTRIBUTIONS_H

#include "utility.h"
#include <stdbool.h>
#include <stdio.h>

typedef struct {
	char key;
	int np;
	real_t * param;
	void * info;
} * Distribution;

void free_Distribution(Distribution dist);
Distribution copy_Distribution(const Distribution dist);
Distribution new_Distribution_from_fp(FILE * fp);
void show_Distribution(FILE * fp, const Distribution dist);
Distribution new_Distribution(const char type, const real_t * param);



int nparameter_distribution(const Distribution dist);
bool validate_parameters(const Distribution dist);
real_t qdistribution(const real_t px, const Distribution dist, const bool tail, const bool logp);
real_t qlogistic(const real_t p, const real_t loc, const real_t sc, const bool tail, const bool logp);

#endif /* LAMBDA_DISTRIBUTIONS_H */
