/*
 * qc_batch.c
 *
 *  Created on: Aug 5, 2011
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "qc_batch.h"

//====================================================================================
//  qc_batch.c
//
//  qc_batch methods for inserting, removing and deleting qc_batch objects
//====================================================================================

char* kmers_string(int index, char* kmer);

//-----------------------------------------------------
// qc_batch_list_init
//
// init to zero the object and initialize the lock
//-----------------------------------------------------

void qc_batch_list_init(qc_batch_list_t *qc_batch_list_p) {
   memset(qc_batch_list_p, 0, sizeof(qc_batch_list_t));
   //qc_batch_list_p->lock = PTHREAD_MUTEX_INITIALIZER;
   pthread_mutex_init(&(qc_batch_list_p->lock), NULL);
}

//-----------------------------------------------------
// qc_batch_list_insert
//
// insert a qc_batch object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

void qc_batch_list_insert(qc_batch_t *qc_batch_p, qc_batch_list_t *qc_batch_list_p) {

  if (qc_batch_list_p==NULL) return;

  pthread_mutex_lock(&qc_batch_list_p->lock);

   //printf("+++++++++++++++++++++++  > qc_batch_list_insert, before, length = %i, qc id = %i\n", qc_batch_list_p->length, qc_batch_p->id);

  if (qc_batch_list_p->first_p==NULL) {
    qc_batch_p->prev_p = NULL;
    qc_batch_p->next_p = NULL;
    qc_batch_list_p->first_p = qc_batch_p;
    qc_batch_list_p->last_p = qc_batch_p;
  } else {
    qc_batch_list_p->last_p->next_p = qc_batch_p;
    qc_batch_p->prev_p = qc_batch_list_p->last_p;
    qc_batch_p->next_p = NULL;
    qc_batch_list_p->last_p = qc_batch_p;
  }
  qc_batch_list_p->length++;

   //printf("+++++++++++++++++++++++  > qc_batch_list_insert, before, length = %i, qc id = %i\n", qc_batch_list_p->length, qc_batch_p->id);

  pthread_mutex_unlock(&qc_batch_list_p->lock);
}

//-----------------------------------------------------
// qc_batch_list_remove
//
// remove the first qc_batch object in the begining
// of the list,
// according to fifo order,
//-----------------------------------------------------

qc_batch_t* qc_batch_list_remove(qc_batch_list_t* qc_batch_list_p) {

  if (qc_batch_list_p==NULL) return NULL;

  pthread_mutex_lock(&qc_batch_list_p->lock);

  //printf("=========================> qc_batch_list_remove, before, length = %i\n", qc_batch_list_p->length);

  // just get the first element, and if is not null update
  // the first-element pointer
  //qc_batch_list_p->length
  qc_batch_t* qc_batch_p = qc_batch_list_p->first_p;
  if (qc_batch_p!=NULL) {
    qc_batch_list_p->first_p = qc_batch_p->next_p;
    qc_batch_list_p->length--;
  }

  //printf("=========================> qc_batch_list_remove, before, after, length = %i, qc id = %i\n", qc_batch_list_p->length, (qc_batch_p==NULL ? -1 : qc_batch_p->id));

  pthread_mutex_unlock(&qc_batch_list_p->lock);

  return qc_batch_p;
}


//-----------------------------------------------------
// qc_batch_list_print
//
// print the content of qc_batch objects in the list
//
//-----------------------------------------------------

void qc_batch_list_print(qc_batch_list_t* list_p) {
  int nt_count;
  int sum_quality;

  pthread_mutex_lock(&list_p->lock);

  printf("Number of items: %i\n", list_p->length);

  qc_batch_t* qc_batch_p = list_p->first_p;

  while(qc_batch_p!=NULL) {
    printf("qc_batch_source_id: %i, qc_batch id: %i\n", qc_batch_p->source_id, qc_batch_p->id);
/*     printf("gpu_quality_result_p:\n");

    for (int j=0; j<200; j++) {
      nt_count = 0;
      sum_quality = 0;

      for (int k=0; k<qc_batch_p->nb_reads; k++) {
	//printf("j: %i, k: %i\n", j, k);
	nt_count += qc_batch_p->gpu_quality_result_p[k].count[j];
	sum_quality += qc_batch_p->gpu_quality_result_p[k].sum_quality[j];
      }

      printf("\tcount[%i]: %i - sum_quality[%i]: %i\n", j, nt_count, j, sum_quality);
    }
*/
    int i, j;
    for(i=0; (i < qc_batch_p->nb_reads && i < 20); i++) {
	    printf("read[%i]: A[%i], C[%i], G[%i], T[%i], N[%i], length: %i\n", i, qc_batch_p->gpu_result_p[i].counters[A], qc_batch_p->gpu_result_p[i].counters[C], qc_batch_p->gpu_result_p[i].counters[G], qc_batch_p->gpu_result_p[i].counters[T], qc_batch_p->gpu_result_p[i].counters[N], qc_batch_p->gpu_result_p[i].counters[TOTAL]);
    }

    printf("kmers from valid reads: \n");
    for(j=0; (j < 20); j++) {
	    printf("qc_batch_p->gpu_kmers_p[%i].kmer: %s, qc_batch_p->gpu_kmers_p[%i].total_count: %i\n", j, qc_batch_p->gpu_kmers_valid_p[j].kmer, j, qc_batch_p->gpu_kmers_valid_p[j].total_count);
    }

    printf("\nkmers from invalid reads: \n");
    for(j=0; (j < 20); j++) {
	    printf("qc_batch_p->gpu_kmers_p[%i].kmer: %s, qc_batch_p->gpu_kmers_p[%i].total_count: %i\n", j, qc_batch_p->gpu_kmers_invalid_p[j].kmer, j, qc_batch_p->gpu_kmers_invalid_p[j].total_count);
    }

    qc_batch_p = qc_batch_p->next_p;
  }

  pthread_mutex_unlock(&list_p->lock);
}

//-----------------------------------------------------
// qc_batch_list_items_free
//
// free list itmes
//
//-----------------------------------------------------

void qc_batch_list_items_free(qc_batch_list_t* list_p) {
  qc_batch_t* item_p;
  char log_message[50];

  while ((item_p = qc_batch_list_remove(list_p)) != NULL) {
    sprintf(log_message, "liberating qc item %i, gpu_result_p = 0x%x\n", item_p->id, item_p->gpu_result_p);
    LOG_DEBUG(log_message);
    qc_batch_free(item_p, 1);
  }
}

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
      int i;
      
      for (i=0; i < KMERS_COMBINATIONS; i++) {
	  qc_kmers_p[i].id = i;
	  qc_kmers_p[i].total_count = 0;
	  memset(qc_kmers_p[i].position_count, 0, MAX_LINE_LENGTH * sizeof(int));
	  //memcpy(&qc_kmers_p[i].kmer, kmers_string(i, kmer), 5 * sizeof(char));
// 	  printf("kmers_string(%i): %s\n", i, kmers_string(i, kmer));
// 	  printf("qc_kmers_p[%i].kmer: %s\n", i, qc_kmers_p[i].kmer);
      }
}

//------------------------------------------------------------------------------------
//  function for kmers string calculation from index
//------------------------------------------------------------------------------------

char* kmers_string(int index, char* kmer) {

	int i, rest, quotient;

	quotient = index;
	kmer[5]='\0';

	for (i=0; i < 5; i++) {
	    rest = quotient % 4;
	    quotient /= 4;

	    switch(rest)
	    {
		case 0 : kmer[4-i] = 'A';
				break;
		case 1 : kmer[4-i] = 'C';
				break;
		case 2 : kmer[4-i] = 'T';
		         break;
		case 3 : kmer[4-i] = 'G';
				break;
	    }

// 	    printf("kmer[%i]: %c\n", i, kmer[4-i]);
// 	    printf("quotient: %i\n", quotient);
// 	    printf("rest: %i\n", rest);
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
