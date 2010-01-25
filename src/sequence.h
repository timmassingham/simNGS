/*
 *  Copyright (C) 2009 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the trim_adapter software for trimming 
 *  adapter sequence from reads.
 *
 *  trim_adapter is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  trim_adapter is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with trim_adapter.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "nuc.h"


struct _sequence {
   char * name, *qname;
   NUC * seq;
   PHREDCHAR * qual;
   uint32_t length;
};
typedef struct _sequence * SEQ;


void free_SEQ( SEQ seq );
SEQ new_SEQ(const uint32_t len, const bool has_qual);
SEQ sequence_from_str( const char * restrict name, const char * restrict seqstr, const char * qname, const char * restrict qualstr );
SEQ copy_SEQ( const SEQ seq);
void show_SEQ(FILE * fp, const SEQ seq);

// Creating
SEQ sequence_from_file  ( FILE * fp);
SEQ sequence_from_fasta ( FILE * fp);
SEQ sequence_from_fastq ( FILE * fp);

//
SEQ reverse_complement_SEQ( const SEQ seq);
#endif

