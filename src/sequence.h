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

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "nuc.h"
#include "utility.h"

typedef struct _cigelt {
	char type;
	int num;
	struct _cigelt * nxt;
} * CIGELT;

typedef struct _ciglist {
	CIGELT start, end;
} CIGLIST;

static CIGLIST null_CIGLIST __attribute__((unused)) = {NULL,NULL};

void show_CIGLIST(FILE * fp, const CIGLIST cigar);
char * string_CIGLIST(const CIGLIST cigar);
CIGLIST copy_CIGLIST(const CIGLIST cigar);
void free_CIGLIST(CIGLIST cigar);
CIGLIST sub_cigar(const CIGLIST cigar, const int len);
CIGLIST reverse_cigar(const CIGLIST cigar);


struct _sequence {
   CIGLIST cigar;
   char * name, *qname;
   ARRAY(NUC) seq;
   ARRAY(PHREDCHAR) qual;
   uint32_t length;
};
typedef struct _sequence * SEQ;


void free_SEQ( SEQ seq );
SEQ new_SEQ(const uint32_t len, const bool has_qual);
SEQ sequence_from_str( const char * restrict name, const char * restrict seqstr, const char * qname, const char * restrict qualstr );
SEQ copy_SEQ( const SEQ seq);
void show_SEQ(FILE * fp, const SEQ seq, CSTRING fmt);

// Creating
SEQ sequence_from_file  ( FILE * fp);
SEQ sequence_from_fasta ( FILE * fp);
SEQ sequence_from_fastq ( FILE * fp);

//
SEQ reverse_complement_SEQ( const SEQ seq, const bool revcigar);
SEQ mutate_SEQ( const SEQ seq, const real_t ins, const real_t del, const real_t mut );
SEQ sub_SEQ( const SEQ seq, const uint32_t loc, const uint32_t len);
#endif

