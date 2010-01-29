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

#ifndef _INTENSITIES_H
#define _INTENSITIES_H

#ifndef _UTILITY_H
#include "utility.h"
#endif

#ifndef _MATRIX_H
#include "matrix.h"
#endif

#ifndef _NUC_H
#include "nuc.h"
#endif

typedef struct {
    uint32_t ncycle,orig_ncycle;
    uint32_t lane,tile;
    bool paired;    // Was read paired?
    MAT cov1,cov2;  // Covariance
    MAT chol1,chol2;       // Cholesky factorisation of covariance
    MAT *invchol1, *invchol2;    // Inverse of cholesky factorisation
    real_t shape, scale; // Parameters for distribution of lambda
    char * label;
} * MODEL;

MODEL new_MODEL(const char * label, const uint32_t lane, const uint32_t tile, const real_t shape, const real_t scale, const MAT cov1, const MAT cov2);
void free_MODEL(MODEL model);
MODEL copy_MODEL(const MODEL model);
void show_MODEL(FILE * fp, const MODEL model);

MODEL trim_MODEL(const uint32_t ncycle, const MODEL model);

MODEL new_MODEL_from_fp( FILE * fp );
MODEL new_MODEL_from_file( const CSTRING filename);

MAT generate_pure_intensities ( const real_t varfact, const real_t lambda, const ARRAY(NUC) seq, const uint32_t ncycle, const MAT chol, MAT ints);
MAT likelihood_cycle_intensities ( const real_t varfact, const real_t mu, const real_t lambda, const MAT ints, const MAT invchol[], MAT like);
void fprint_intensities(FILE * fp, const char * prefix, const MAT ints, const bool last);
NUC * call_by_maximum_likelihood(const MAT likelihood, NUC * calls);
uint32_t number_inpure_cycles( const MAT intensities, const real_t threshold, const uint32_t ncycle);

#endif
