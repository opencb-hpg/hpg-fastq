/*
 * Copyright (c) 2012 Victor Requena (BULL)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
 *
 * This file is part of hpg-fastq.
 *
 * hpg-fastq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hpg-fastq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-fastq. If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef QC_BATCH_H
#define QC_BATCH_H

#include "commons.h"
#include "fastq/fastq_batch.h"

#define A 				1	//Numer of A's: ltrim (10), rtrim (10), total (12)
#define C 				3	//Numer of C's: ltrim (10), rtrim (10), total (12)
#define G 				7	//Numer of G's: ltrim (10), rtrim (10), total (12)
#define N 				6	//Numer of N's: ltrim (10), rtrim (10), total (12)
#define T 				4	//Numer of T's: ltrim (10), rtrim (10), total (12)
#define TOTAL 				2
#define MEAN_MEDIAN_QUALITY 		0	//Median (8) / Mean (8) in the defined quality window
#define SIDE_NTS_QUALITY 		5	//Mean quality of the left and right sides for trim and filter: ltrim (8) / lfilter (8) / rfilter (8) / ltrim (8)
#define TRIM_QUALITY 			8	//Mean quality of the trimmed read: ltrim (8) / rtrim (8) / trim (8) / total (8)
#define OUT_OF_QUALITY			9	//Number of nucleotides out of quality thresholds: left (10), right (10), no cut (12)

#define GET_0_8(x)			(x & 255)
#define GET_8_16(x)			((x >> 8) & 255)
#define GET_16_24(x)			((x >> 16) & 255)
#define GET_24_32(x)			((x >> 24) & 255)
#define GET_0_12(x)			(x & 4095)
#define GET_12_22(x)			((x >> 12) & 1023)
#define GET_22_32(x)			((x >> 22) & 1023)

#define COUNTERS_SIZE 			10
#define INT_COUNTERS_SIZE_IN_MEMORY 	40
#define MAX_LINE_LENGTH 		1024
#define KMERS_COMBINATIONS		1024

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Structure for QC counters of a read 
* 
* Vector containing QC values for a given read (#A, #C, #G, #T, #N, length, ...)
*/
typedef struct qc_read {
    int counters[COUNTERS_SIZE];		/**< Vector containing A, C, G, T, N count, mean&median quality */
} qc_read_t;

/**
* @brief Structure for QC quality by nt position
* 
* Vector containing accumulated values of quality by nt position
*/
typedef struct qc_quality {
    int sum_quality[MAX_LINE_LENGTH];	/**< Vector containing the quality sum by read position */
    int count[MAX_LINE_LENGTH];		/**< Vector containing the count of nts by read position */
} qc_quality_t;

/**
* @brief Structure for kmers statistics
* 
* Structure containing kmers statistics for a given kmer: id, kmer string, kmer count, kmer count by position
*/
typedef struct qc_kmers {
    int id;						/**< Kmer id */
    char kmer[5];					/**< String representation of the kmer (e.g.: AATCG) */
    unsigned int total_count;				/**< Total count of the kmer */
    unsigned int position_count[MAX_LINE_LENGTH];	/**< Count of the kmer by nt position */
} qc_kmers_t;

/**
* @brief Batch of QC information
* 
* Structure containing quality control information. A qc_batch contains information for a fastq_batch.
*/
typedef struct qc_batch {
    int id;					/**< Batch id */
    int source_id;				/**< Source id (pair 1 or 2) */
    int nb_reads;				/**< Number of reads in the batch */
    fastq_batch_t* read_p;			/**< Pointer to batch of fastq reads */
    qc_read_t* gpu_result_p;			/**< Pointer to qc_read_t results */
    qc_quality_t* gpu_quality_result_p;		/**< Pointer to qc_quality_t results */
    int* gpu_nt_type_valid_counter_p;		/**< Pointer to a nt position counter of valid reads */
    int* gpu_nt_type_invalid_counter_p;		/**< Pointer to a nt position counter of invalid reads */
    qc_kmers_t* gpu_kmers_valid_p;		/**< Pointer to qc_kmers structure of valid reads */
    qc_kmers_t* gpu_kmers_invalid_p;		/**< Pointer to qc_kmers structure of invalid reads */
    struct qc_batch* prev_p;			/**< Pointer to the previous qc_batch */
    struct qc_batch* next_p;			/**< Pointer to the next qc_batch */
} qc_batch_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Frees a qc_batch structure
*  @param[in,out] qc_batch_p pointer to the qc_batch to be freed
*  @param all indicates if reads data must be freed
*  @return void
*  
*  This function frees a qc_batch structure with all = 1 frees read data, with all = 0 preserves read data
*/
void qc_batch_free(qc_batch_t* qc_batch_p, int all);

/**
*  @brief Inits a qc_kmer structure
*  @param[in,out] qc_kmers_p pointer to the qc_kmer structure
*  @return void
*  
*  This function initializes a qc_kmer structure
*/
void qc_kmers_init(qc_kmers_t* qc_kmers_p);

/**
*  @brief Returns the string representation of a kmer given its id
*  @param index pointer to the qc_batch to be freed
*  @param kmer pointer to the kmer string
*  @return char* pointer to the filled kmer string
*  
*  This function returns the string representation of a kmer provided its id
*  For example, id = 0 -> kmer = "AAAAA", id = 1023 -> kmer = "TTTTT"
*/
char* kmers_string(int index, char* kmer);

/**
*  @brief Auxiliar function for sorting kmer structures
*  @param k1 kmer structure 1
*  @param k2 kmer structure 2
*  @return int
*  
*  Returns > 0 if kmer2 total count is greater than kmer1, 0 if both total count of the kmers are equal, and < 0 if 
*  kmer1 total count is greater than kmer2
*/
int kmers_sort(const void* k1, const void* k2);

#endif /* QC_BATCH_H */
