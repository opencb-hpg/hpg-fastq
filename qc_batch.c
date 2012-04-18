
#include <stdlib.h>
#include <pthread.h>
#include <string.h>

#include "qc_batch.h"

//====================================================================================
//  qc_batch.c
//
//  qc_batch methods for inserting, removing and deleting qc_batch objects
//====================================================================================

char* kmers_string(int index, char* kmer);

//-----------------------------------------------------
// qc_batch_free
//
// free memory for this structure,
// the input parameter 'all' indicates free memory
// associated to the pointers in his structure
//-----------------------------------------------------

void qc_batch_free(qc_batch_t* qc_batch_p, int all) {
  if (qc_batch_p==NULL) return;

  if (all) {
    if (qc_batch_p->read_p!=NULL) fastq_batch_free(qc_batch_p->read_p); //free(qc_batch_p->read_p);
    if (qc_batch_p->gpu_result_p!=NULL) free(qc_batch_p->gpu_result_p);
    if (qc_batch_p->gpu_nt_type_valid_counter_p!=NULL) free(qc_batch_p->gpu_nt_type_valid_counter_p);
    if (qc_batch_p->gpu_nt_type_invalid_counter_p!=NULL) free(qc_batch_p->gpu_nt_type_invalid_counter_p);
    if (qc_batch_p->gpu_kmers_valid_p!=NULL) free(qc_batch_p->gpu_kmers_valid_p);
    if (qc_batch_p->gpu_kmers_invalid_p!=NULL) free(qc_batch_p->gpu_kmers_invalid_p);
  }
  free(qc_batch_p);
}

//-----------------------------------------------------
// qc_kmers_init
//
// initialize values of kmers indexes and strings
//-----------------------------------------

void qc_kmers_init(qc_kmers_t* qc_kmers_p) {  
  char kmer[6];
  
  for (int i=0; i<KMERS_COMBINATIONS; i++) {
    qc_kmers_p[i].id = i;
    qc_kmers_p[i].total_count = 0;
    memset(qc_kmers_p[i].position_count, 0, MAX_LINE_LENGTH * sizeof(int));
  }
}

//------------------------------------------------------------------------------------
//  function for kmers string calculation from index
//------------------------------------------------------------------------------------

char* kmers_string(int index, char* kmer) {
	int rest, quotient;

	quotient = index;
	kmer[5]='\0';

	for (int i=0; i<5; i++) {
	    rest = quotient % 4;
	    quotient /= 4;

	    switch(rest) {
		case 0: 
		    kmer[4-i] = 'A';
		    break;
		case 1: 
		    kmer[4-i] = 'C';
		    break;
		case 2: 
		    kmer[4-i] = 'T';
		    break;
		case 3: 
		    kmer[4-i] = 'G';
		    break;
	    }
	}

	return kmer;
}

//------------------------------------------------------------------------------------
//  function for kmer count sorting
//------------------------------------------------------------------------------------

int kmers_sort(const void* k1, const void* k2) {
      qc_kmers_t* kmer1_p = (qc_kmers_t*) k1;
      qc_kmers_t* kmer2_p = (qc_kmers_t*) k2;

      return (kmer2_p->total_count - kmer1_p->total_count);
}