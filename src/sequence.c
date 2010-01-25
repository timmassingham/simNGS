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


void free_SEQ ( SEQ seq ){
   validate(NULL!=seq,);
   if ( NULL!=seq->name ){ safe_free(seq->name); }
   if ( NULL!=seq->qname){ safe_free(seq->qname);}
   if ( NULL!=seq->seq  ){ safe_free(seq->seq); }
   if ( NULL!=seq->qual ){ safe_free(seq->qual); }
   safe_free(seq);
}


SEQ new_SEQ (const uint32_t len, const bool has_qual){
   SEQ seq = calloc(1,sizeof(struct _sequence));
   if(NULL==seq){return NULL;}
   seq->name = NULL;
   seq->qname = NULL;
   seq->seq = calloc(len,sizeof(NUC));
   if(NULL==seq->seq){ goto cleanup; }
   if(has_qual){
      seq->qual = calloc(len,sizeof(PHREDCHAR));
      if(NULL==seq->qual){ goto cleanup; }
   }

   seq->length = len;

   return seq;

cleanup:
   free_SEQ(seq);
   return NULL;
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

   for ( uint32_t i=0 ; i<length ; i++){
      seq->seq[i] = nuc_from_char(seqstr[i]);
   }
   if ( NULL!=qualstr ){
      for ( uint32_t i=0 ; i<length ; i++){
         seq->qual[i] = phredchar_from_char(qualstr[i]);
      }
   }
   
   return seq;
}

SEQ copy_SEQ ( const SEQ seq){
   if(NULL==seq){ return NULL;}

   SEQ newseq = calloc(1,sizeof(struct _sequence));
   if(NULL==newseq){return NULL;}
   // Sequence
   newseq->seq = calloc(seq->length,sizeof(NUC));
   if(NULL==newseq->seq){ goto cleanup;}
   memcpy(newseq->seq,seq->seq,seq->length*sizeof(NUC));
   // Qualities
   if(NULL!=seq->qual){
       newseq->qual = calloc(seq->length,sizeof(PHREDCHAR));
       if(NULL==newseq->qual){ goto cleanup;}
       memcpy(newseq->qual,seq->qual,seq->length*sizeof(PHREDCHAR));
   }
   
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

   (NULL==seq->qual)? fputc('>',fp) : fputc('@',fp);
   if(NULL!=seq->name){ fputs(seq->name,fp); }
   fputc('\n',fp);
   for ( uint32_t i=0 ; i<seq->length ; i++){
       show_NUC(fp,seq->seq[i]);
   }
   fputc('\n',fp);
   if(NULL!=seq->qual){
      fputc('+',fp);
      if(NULL!=seq->qname){ fputs(seq->qname,fp); }
      fputc('\n',fp);
      fwrite(seq->qual,sizeof(char),seq->length,fp);
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

/*  Utility function. Skip until find character
 * Returns character on success or EOF if character
 * not found before end of file
 */
int skipUntilChar ( FILE * fp, const char c){
   assert(NULL!=fp);

   int n;
   do {
      n = fgetc(fp);
   } while (EOF!=n && c!=n);
   return n;
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
    NUC * rcnuc = NULL;
    SEQ newseq = NULL;
    
    validate(NULL!=seq,NULL);
    newseq = copy_SEQ(seq);
    validate(NULL!=newseq,NULL);
    rcnuc = reverse_complement(seq->seq,seq->length);
    if(NULL==rcnuc){ goto cleanup;}
    safe_free(newseq->seq);
    newseq->seq = rcnuc;
    
    return newseq;
    
cleanup:
    safe_free(rcnuc);
    free_SEQ(newseq);
    return NULL;
}

#ifdef TEST
int main(int argc, char * argv[]){
    FILE * fp = (argc==1)?stdin:fopen(argv[1],"r");
    
    SEQ seq=NULL;
    while( NULL!=(seq=sequence_from_file(fp)) ){
        show_SEQ(stdout,seq);
        SEQ rcseq = reverse_complement_SEQ(seq);
        show_SEQ(stdout,rcseq);
        free_SEQ(rcseq);
        free_SEQ(seq);
    }
    return EXIT_SUCCESS;
}
#endif
