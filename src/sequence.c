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

#include "nuc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <err.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include "utility.h"
#include "sequence.h"
#include "mystring.h"
#include "random.h"

void free_CIGLIST(CIGLIST cigar){
        CIGELT elt = cigar.start;
        while(NULL!=elt){
                CIGELT nxtelt = elt->nxt;
                free(elt);
                elt = nxtelt;
        }
}

// Drops deletions from beginning of cigar string
CIGLIST compact_CIGLIST(CIGLIST cigar){
	CIGELT elt = cigar.start;
	if(NULL==elt){ return cigar; }

	if('D'==elt->type){
		CIGELT nxtelt = elt->nxt;
		free(elt);
		if(cigar.start==cigar.end){
			cigar.end = nxtelt;
		}
		cigar.start = nxtelt;
	}
	return cigar;
}

CIGLIST pushStart_CIGLIST(CIGLIST cigar, const char type, const int num){
	CIGELT elt = malloc(sizeof(*elt));
	if(NULL==elt){ 
		free_CIGLIST(cigar);
		return null_CIGLIST;
	}

	elt->nxt = cigar.start;
	elt->type = type;
	elt->num = num;

	cigar.start = elt;
	if(NULL==cigar.end){ cigar.end = elt; }
	return cigar;
}

CIGLIST pushEnd_CIGLIST(CIGLIST cigar, const char type, const int num){
        CIGELT elt = malloc(sizeof(*elt));
        if(NULL==elt){ 
		free_CIGLIST(cigar);
		return null_CIGLIST;
	}

        elt->nxt = NULL;
        elt->type = type;
        elt->num = num;

	if(NULL==cigar.end){
		cigar.start = elt;
		cigar.end = elt;
	} else {
		cigar.end->nxt = elt;
		cigar.end = elt;
	}
	return cigar;
}

void show_CIGLIST(FILE * fp, const CIGLIST cigar){
	if(NULL==fp){ return; }
	CIGELT elt = cigar.start;
        while(NULL!=elt){
		fprintf(fp,"%d%c",elt->num,elt->type);
		elt = elt->nxt;
	}
}

CIGLIST copy_CIGLIST(const CIGLIST cigar){
	CIGLIST newcigar = {NULL,NULL};
	CIGELT elt = cigar.start;

	while(NULL!=elt){
		newcigar = pushEnd_CIGLIST(newcigar,elt->type,elt->num);
		elt = elt->nxt;
	}
	return newcigar;
}

CIGLIST reverse_cigar(const CIGLIST cigar){
	CIGLIST newcig = {NULL,NULL};
	CIGELT elt = cigar.start;
	while(NULL!=elt){
		newcig = pushStart_CIGLIST(newcig,elt->type,elt->num);
		elt = elt->nxt;
	}
	return newcig;
}

CIGLIST sub_cigar(const CIGLIST cigar, const int len){
	CIGLIST newcig = {NULL,NULL};
	CIGELT elt = cigar.start;
	int tot=0;
	while(NULL!=elt && tot<len){
		newcig = pushEnd_CIGLIST(newcig,elt->type,elt->num);
		if('D'!=elt->type){tot += elt->num;}
		elt = elt->nxt;
	}

	if(tot<len){
		// Invalid shortening, pad with N's
		newcig = pushEnd_CIGLIST(newcig,'N',len-tot);
	} else {
		if(NULL!=newcig.end){newcig.end->num -= (tot-len);}
	}

	return newcig;
}



bool hasQual( const SEQ seq){
    validate(NULL!=seq,false);
    return (NULL!=seq->qual.elt)?true:false;
}


void free_SEQ ( SEQ seq ){
   validate(NULL!=seq,);
   if ( NULL!=seq->name ){ safe_free(seq->name); }
   if ( NULL!=seq->qname){ safe_free(seq->qname);}
   if ( NULL!=seq->seq.elt  ){ free_ARRAY(NUC)(seq->seq); }
   if ( NULL!=seq->qual.elt ){ free_ARRAY(PHREDCHAR)(seq->qual); }
   free_CIGLIST(seq->cigar);
   safe_free(seq);
}


SEQ new_SEQ (const uint32_t len, const bool has_qual){
   SEQ seq = calloc(1,sizeof(struct _sequence));
   if(NULL==seq){return NULL;}
   seq->name = NULL;
   seq->qname = NULL;
   seq->cigar = null_CIGLIST;
   seq->seq = new_ARRAY(NUC)(len);
   if(NULL==seq->seq.elt){ goto cleanup; }
   if(has_qual){
      seq->qual = new_ARRAY(PHREDCHAR)(len);
      if(NULL==seq->qual.elt){ goto cleanup; }
   }

   seq->length = len;

   return seq;

cleanup:
   free_SEQ(seq);
   return NULL;
}

SEQ resize_SEQ( SEQ seq, const uint32_t newlen){
    seq->seq = resize_ARRAY(NUC)(seq->seq,newlen);
    if(0==seq->seq.nelt){ free_SEQ(seq); return NULL;}
    for ( uint32_t i=seq->length ; i<newlen ; i++){
        seq->seq.elt[i] = NUC_AMBIG;
    }
    CIGLIST newcig = sub_cigar(seq->cigar,newlen);
    free_CIGLIST(seq->cigar);
    seq->cigar = newcig;
    if(hasQual(seq)){
        seq->qual = resize_ARRAY(PHREDCHAR)(seq->qual,newlen);
        if(0==seq->qual.nelt){ free_SEQ(seq); return NULL;}
        for ( uint32_t i=seq->length ; i<newlen ; i++){
            seq->qual.elt[i] = MIN_PHRED;
        }
    }
    seq->length = newlen;
    return seq;
}

SEQ sequence_from_str ( const char * restrict name, const char * restrict seqstr, const char * restrict qname, const char * restrict qualstr ){
   if(NULL==seqstr){ return NULL;}
   const uint32_t length = strlen(seqstr);
   if ( NULL!=qualstr && strlen(qualstr)!=length){ return NULL;}
   SEQ seq = new_SEQ(length,(NULL!=qualstr)?true:false);
   if(NULL==seq){ return NULL;}
   if(NULL!=name){
      seq->name = calloc(strlen(name)+1,sizeof(char));
      if(NULL==seq->name){ free_SEQ(seq); return NULL;}
      strcpy(seq->name,name);
   }
   if(NULL!=qname){
      seq->qname = calloc(strlen(qname)+1,sizeof(char));
      if(NULL==seq->qname){ free_SEQ(seq); return NULL;}
      strcpy(seq->qname,qname);
   }
   free_CIGLIST(seq->cigar);
   seq->cigar = pushStart_CIGLIST(seq->cigar,'M',length);

   for ( uint32_t i=0 ; i<length ; i++){
      seq->seq.elt[i] = nuc_from_char(seqstr[i]);
   }
   if ( NULL!=qualstr ){
      for ( uint32_t i=0 ; i<length ; i++){
         seq->qual.elt[i] = phredchar_from_char(qualstr[i]);
      }
   }
   
   return seq;
}

SEQ copy_SEQ ( const SEQ seq){
   if(NULL==seq){ return NULL;}

   SEQ newseq = new_SEQ(seq->length,hasQual(seq)?true:false);
   if(NULL==newseq){return NULL;}
   // Sequence
   memcpy(newseq->seq.elt,seq->seq.elt,seq->length*sizeof(NUC));
   // Qualities
   if(hasQual(seq)){
       memcpy(newseq->qual.elt,seq->qual.elt,seq->length*sizeof(PHREDCHAR));
   }
   // cigar
   free_CIGLIST(newseq->cigar);
   newseq->cigar = copy_CIGLIST(seq->cigar);
   
   newseq->length = seq->length;
   if(NULL!=seq->name){
      newseq->name = calloc(strlen(seq->name)+1,sizeof(char));
      if ( NULL==newseq->name){ goto cleanup;}
      strcpy(newseq->name,seq->name);
   }
   if(NULL!=seq->qname){
      newseq->qname = calloc(strlen(seq->qname)+1,sizeof(char));
      if ( NULL==newseq->qname){ goto cleanup;}
      strcpy(newseq->qname,seq->qname);
   }

   return newseq;
   
cleanup:
    free_SEQ(newseq);
    return NULL;
}


void show_SEQ( FILE * fp, const SEQ seq){
   validate(NULL!=fp,);
   validate(NULL!=seq,);

   (!hasQual(seq)) ? fputc('>',fp) : fputc('@',fp);
   if(NULL!=seq->name){ fputs(seq->name,fp); fputc(' ',fp);}
   show_CIGLIST(fp,seq->cigar);
   fputc('\n',fp);
   for ( uint32_t i=0 ; i<seq->length ; i++){
       show_NUC(fp,seq->seq.elt[i]);
   }
   fputc('\n',fp);
   if(hasQual(seq)){
      fputc('+',fp);
      if(NULL!=seq->qname){ fputs(seq->qname,fp); }
      fputc('\n',fp);
      fwrite(seq->qual.elt,sizeof(char),seq->length,fp);
      fputc('\n',fp);
   }
}




/*
 * Functions for reading in sequence
 */

char * read_until(FILE * fp, const char target, const bool skipSpace){
   validate(NULL!=fp,NULL);
   Mystring str = new_mystring(80);
   validate(NULL!=str,NULL);
   char c = fgetc(fp);
   if ( EOF==c ){
      free_mystring(str);
      return NULL;
   }
   while( c!=EOF && c!=target ){
      if( !skipSpace || !isspace(c) ){ append_char_to_mystring (c,str); }
      c = fgetc(fp);
   } 
   char * cstr = cstring_of_mystring(str);
   free_mystring(str);
   return cstr;
}


SEQ sequence_from_fasta ( FILE * fp){
   char *name=NULL, *cseq=NULL;
   SEQ seq = NULL;

   validate(NULL!=fp,NULL);

   // Must start with a '>'
   //int a = fgetc(fp);
   if ( '>' != fgetc(fp) ){ goto sequence_from_fasta_error; }
   // Read name (for now)
   name = read_until(fp,'\n',false);
   if (NULL==name){ goto sequence_from_fasta_error; }
   // Read sequence. assuming all on one line
   cseq = read_until(fp,'>',true);
   ungetc('>',fp);
   if (NULL==cseq){goto sequence_from_fasta_error; }

   seq = sequence_from_str(name,cseq,NULL,NULL);

sequence_from_fasta_error:

   if(NULL!=name){free(name);}
   if(NULL!=cseq){free(cseq);}

   return seq;
}

SEQ sequence_from_fastq ( FILE * fp){
   char *name=NULL, *cseq=NULL, *qname=NULL, *qual=NULL;
   SEQ seq=NULL;

   validate(NULL!=fp,NULL);

   // Must start with a '>'
   if ( '@' != fgetc(fp) ){ goto sequence_from_fastq_error;}
   // Read name (for now)
   name = read_until(fp,'\n',false);
   if (NULL==name){ goto sequence_from_fastq_error; }
   // Read sequence
   cseq = read_until(fp,'+',true);
   if (NULL==cseq){goto sequence_from_fastq_error; }
   // Read quality name
   qname = read_until(fp,'\n',false);
   if (NULL==qname){ goto sequence_from_fastq_error; }
   // Read qualities. assuming all on one line
   qual = read_until(fp,'@',true);
   ungetc('@',fp);
   if (NULL==qual){goto sequence_from_fastq_error; }

   seq = sequence_from_str(name,cseq,qname,qual);

sequence_from_fastq_error:
   if(NULL!=name){free(name);}
   if(NULL!=cseq){free(cseq);}
   if(NULL!=qname){free(qname);}
   if(NULL!=qual){free(qual);}
   
   return seq;
}

SEQ sequence_from_file ( FILE * fp){
   validate(NULL!=fp,NULL);
   int c = fgetc(fp);
   ungetc(c,fp);
   switch(c){
      case '>': return sequence_from_fasta(fp);
      case '@': return sequence_from_fastq(fp);
   }
   return NULL;
}

SEQ reverse_complement_SEQ( const SEQ seq){
    ARRAY(NUC) rcnuc = null_ARRAY(NUC);
    ARRAY(PHREDCHAR) revqual = null_ARRAY(PHREDCHAR);
    SEQ newseq = NULL;
    
    validate(NULL!=seq,NULL);
    newseq = copy_SEQ(seq);
    validate(NULL!=newseq,NULL);

    rcnuc = reverse_complement(seq->seq);
    if(NULL==rcnuc.elt){ goto cleanup;}
    free_ARRAY(NUC)(newseq->seq);
    newseq->seq = rcnuc;

    if(NULL!=seq->qual.elt){
        revqual = reverse_quality(seq->qual);
        if(NULL==revqual.elt){ goto cleanup; }
        free_ARRAY(PHREDCHAR)(newseq->qual);
        newseq->qual = revqual;
    }

    free_CIGLIST(newseq->cigar);
    newseq->cigar = reverse_cigar(seq->cigar);
    
    return newseq;
    
cleanup:
    free_ARRAY(NUC)(rcnuc);
    free_ARRAY(PHREDCHAR)(revqual);
    free_SEQ(newseq);
    return NULL;
}


SEQ mutate_SEQ ( const SEQ seq, const real_t ins, const real_t del, const real_t mut ){
    validate(NULL!=seq,NULL);
    validate(isprob(ins),NULL);
    validate(isprob(del),NULL);
    validate(isprob(mut),NULL);
    const real_t probs[4] = { ins, del, mut, 1.-ins-del-mut };
    
    SEQ mutseq = copy_SEQ(seq);
    uint32_t scount=0, mcount=0;
    char cigType = 0; int cigNum = 0;
    CIGLIST cigar = null_CIGLIST;
    while(scount<seq->length){
        if(mcount==mutseq->length){ // Enlarge mutated sequence
            resize_SEQ(mutseq,mutseq->length*2);
        }
        uint32_t i = rchoose(probs,4);
        switch(i){
           case 0: // Insertion of base
	       if('I'==cigType){ cigNum++;}
	       else {
		       if(0!=cigType){cigar = pushEnd_CIGLIST(cigar,cigType,cigNum);}
		       cigType = 'I';
		       cigNum = 1;
	       }
               mutseq->seq.elt[mcount] = random_NUC();
               mcount++;
               break;
           case 1: // Deletion
	       if('D'==cigType){ cigNum++;}
	       else {
		      if(0!=cigType){cigar = pushEnd_CIGLIST(cigar,cigType,cigNum);}
		      cigType = 'D';
		      cigNum = 1;
	       } 
               scount++;
               break;
           case 2: // Mutation
	       if('S'==cigType){ cigNum++;}
	       else {
                      if(0!=cigType){cigar = pushEnd_CIGLIST(cigar,cigType,cigNum);}
                      cigType = 'S';
                      cigNum = 1;
               }
               mutseq->seq.elt[mcount] = random_other_NUC(seq->seq.elt[scount]);
               mcount++; scount++;
               break;
           case 3: // Match
	       if('M'==cigType){ cigNum++;}
	       else {
                      if(0!=cigType){cigar = pushEnd_CIGLIST(cigar,cigType,cigNum);}
                      cigType = 'M';
                      cigNum = 1;
               }
               mutseq->seq.elt[mcount] = seq->seq.elt[scount];
               mcount++; scount++;
               break;
           default:
               errx(EXIT_FAILURE,"Invalid case '%d' generated in %s (%s:%u)",i,__func__,__FILE__,__LINE__);
        }
    }
    // Push any remaining edits into cigar string, except deletions.
    if(0!=cigType && 'D'!=cigType){cigar = pushEnd_CIGLIST(cigar,cigType,cigNum);}
    cigar = compact_CIGLIST(cigar);
    mutseq = resize_SEQ(mutseq,mcount);
    free_CIGLIST(mutseq->cigar);
    mutseq->cigar = cigar;
    return mutseq;
}

SEQ sub_SEQ( const SEQ seq, const uint32_t loc, const uint32_t len){
    validate(NULL!=seq,NULL);
    
    SEQ subseq = new_SEQ(len,hasQual(seq));
    validate(NULL!=subseq,NULL);
    
    const uint32_t seqend = (loc+len<seq->length)?(len):(seq->length-loc);
    // Name
    subseq->name = copy_CSTRING(seq->name);
    // Sequence
    for ( uint32_t i=0 ; i<seqend ; i++){
        subseq->seq.elt[i] = seq->seq.elt[loc+i];
    }
    for( uint32_t i=seqend ; i<len ; i++){
        subseq->seq.elt[i] = NUC_AMBIG;
    }
    // Qualities
    if(hasQual(seq)){
        for ( uint32_t i=0 ; i<seqend ; i++){
            subseq->qual.elt[i] = seq->qual.elt[loc+i];
        }
        for( uint32_t i=seqend ; i<len ; i++){
            subseq->qual.elt[i] = MIN_PHRED;
        }
    }
    // Cigar string
    free_CIGLIST(subseq->cigar);
    subseq->cigar = sub_cigar(seq->cigar,len);

    return subseq;
}

#ifdef TEST
int main(int argc, char * argv[]){
    FILE * fp = (argc==1)?stdin:fopen(argv[1],"r");
    init_gen_rand(12351);
    
    SEQ seq=NULL;
    while( NULL!=(seq=sequence_from_file(fp)) ){
        show_SEQ(stdout,seq);
	SEQ mutseq = mutate_SEQ(seq,0.02,0.02,0.02);
	show_SEQ(stdout,mutseq);
        SEQ rcseq = reverse_complement_SEQ(mutseq);
        show_SEQ(stdout,rcseq);
        free_SEQ(rcseq);
	free_SEQ(mutseq);
        free_SEQ(seq);
    }
    return EXIT_SUCCESS;
}
#endif
