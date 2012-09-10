
#ifndef PREPRO_BATCH_H
#define PREPRO_BATCH_H

#include <pthread.h>

#include "commons.h"
#include "fastq/fastq_batch.h"
#include "fastq/fastq_file.h"

#define VALID_READ		'V'
#define INVALID_READ	 	'I'
#define TRIM_READ	        'T'
#define LTRIM_READ	        'L'
#define RTRIM_READ	        'R'
#define UNKNOWN_STATUS_READ     'U'

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Batch with status information about reads
* 
* Structure containing the status of the reads for filtering and preprocessing purposes
*/
typedef struct status_batch {
    int id;				/**< Id of the batch. */
    int source_id;			/**< Source id (pairend 1 or 2) of the batch. */
    int num_reads;			/**< Number of reads of the batch. */
    char* read_status;			/**< Read status (valid, invalid). */
    fastq_batch_t* read_p;		/**< Pointer to fastq reads. */
    struct status_batch* prev_p;	/**< Pointer to the previous batch. */
    struct status_batch* next_p;	/**< Pointer to the next batch. */
} status_batch_t;

/**
* @brief List of status batch
* 
* List to contain and link status batches
*/
typedef struct status_batch_list {
    int length;					/**< Total length of the list. */
    int length_by_source_id[MAX_NUM_PRODUCERS];	/**< Length of the list by source id. */
    pthread_mutex_t lock;			/**< Lock.  */
    struct status_batch* first_p;		/**< Pointer to the first batch in the list. */
    struct status_batch* last_p;		/**< Pointer to the last batch in the list. */
} status_batch_list_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Init a status batch list
*  @param[in,out] status_batch_list_p pointer to the list to initialize
*  @return void
*  
*  Performs status batch list initialization
*/
void status_batch_list_init(status_batch_list_t* status_batch_list_p);

/**
*  @brief Inserts status batch in a list
*  @param status_batch_p pointer to the status batch to insert in the list
*  @param[in,out] status_batch_list_p pointer to the list in where batch will be inserted
*  @return void
*  
*  Inserts a status batch in a status batch list
*/
void status_batch_list_insert(status_batch_t* status_batch_p, status_batch_list_t* status_batch_list_p);

/**
*  @brief Removes status batch from a list
*  @param[in,out] status_batch_list_p pointer to the list from where batch will be removed
*  @return status_batch_t pointer to the removed batch
*  
*  Removes the first batch from the list and returns it
*/
status_batch_t* status_batch_list_remove(status_batch_list_t* status_batch_list_p);

/**
*  @brief Removes status batch from a list with the specified source_id and batch_id
*  @param[in,out] status_batch_list_p pointer to the list from where batch will be removed
*  @param source_id source id (pairend 1 or 2) of the batch to remove
*  @param batch_id  id of the batch to remove
*  @return status_batch_t pointer to the removed batch
*  
*  Removes the specified batch from the list and returns it
*/
status_batch_t* status_batch_list_get_next_by_id(status_batch_list_t* status_batch_list_p, int source_id, int batch_id);

/**
*  @brief Prints a status batch list in the console
*  @param status_batch_list_p pointer to the status batch list
*  @return void
*  
*  Prints the content of each status batch from the list in the console
*/
void status_batch_list_print(status_batch_list_t* status_batch_list_p);

/**
*  @brief Frees status batch list items
*  @param[in,out] status_batch_list_p pointer to the status batch list
*  @return void
*  
*  Free status batch list items
*/
void status_batch_list_items_free(status_batch_list_t* list_p);

/**
*  @brief Frees status batch 
*  @param[in,out] status_batch_p pointer to the status batch
*  @param all flag for indicate read data free
*  @return void
*  
*  Free status batch, with all = 1 frees read data, with all = 0 preserves read data
*/
void status_batch_free(status_batch_t* status_batch_p, int all);

/**
*  @brief Calculates the mean quality of the last nucleotides of a read
*  @param quality qualities of the nucleotides (all batch reads qualities)
*  @param data_indices indices of data&quality to index the quality positions in the vector
*  @param position position of the read in the batch
*  @param last_nts number of nucleotides to compute for mean quality
*  @return int mean quality of the last specified nucleotides
*  
*  Calculate the specified last nucleotides mean quality for a given read within
*  the vector
*/
int calculate_last_nt_mean_quality(char* quality, int* data_indices, int position, int last_nts);

#endif /* PREPRO_BATCH_H */
