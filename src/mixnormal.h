#ifndef MIXNORMAL_H
#define MIXNORMAL_H

typedef struct {
        unsigned int nmix;
	real_t * prob;
        real_t * mean;
        real_t * sd;
} * NormMixParam;


NormMixParam fit_mixnormal(const real_t * x, const unsigned int n, const unsigned int nmix, const unsigned niter);
NormMixParam new_NormMixParam( const unsigned int nmix);
void free_NormMixParam(NormMixParam param);
void show_NormMixParam( FILE * fp, const NormMixParam param);
NormMixParam copy_NormMixParam(const NormMixParam param);


real_t dmixnorm(const real_t x, const NormMixParam param, const bool logd);
real_t rmixnorm(const NormMixParam param);
real_t pmixnorm(const real_t x, const NormMixParam param, const bool tail, const bool logp);
real_t qmixnorm(const real_t p, const NormMixParam param, const bool tail, const bool logp);

#endif /* MIXNORMAL_H */

