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


#ifndef PREPRO_H
#define PREPRO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>

#include "commons.h"
#include "fastq/fastq_batch_list.h"
#include "fastq/fastq_file.h"
#include "list.h"
#include "log.h"
#include "file_utils.h"
#include "prepro_commons.h"
#include "system_utils.h"

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief QC calculator server parameters
* 
* Structure containing parameters to pass to the QC calculator server
*/
// typedef struct qc_calc_server_input {
//     int num_gpu_devices;		/**< Number of GPU devices. */
//     int cpu_num_threads;		/**< Number of CPU threads. */
//     int gpu_device_id[256];		/**< GPU device identifiers to use. */
//     int num_blocks;			/**< Number of GPU blocks to launch. */
//     int num_threads;			/**< Number of GPU threads by warp. */
//     int min_quality;			/**< Minimum accepted quality. */
//     int max_quality;			/**< Maximum accepted quality. */
//     int begin_quality_nt;		/**< First nt in the read to compute mean and median quality. */
//     int end_quality_nt;			/**< Last nt in the read to compute mean and median quality. */
//     int rtrim_nts;			/**< Number of last nts to compute statistics for preprocessing. */
//     int ltrim_nts;			/**< Number of first nts to compute statistics for preprocessing. */
//     int rfilter_nts;			/**< Number of last nts to compute statistics for filtering. */
//     int lfilter_nts;			/**< Number of first nts to compute statistics for filtering. */
//     int kmers_on;			/**< Flag for kmers calculation. */
//     fastq_batch_list_t* batch_list_p;	/**< Pointer to the fastq batch list. */
// } qc_calc_server_input_t;

/**
* @brief Results server parameters
* 
* Structure containing parameters to pass to the results server
*/
// typedef struct results_server_input {
//     int num_blocks;			/**< Number of GPU blocks launched. */
//     int num_threads;			/**< Number of GPU threads launched by warp. */
//     int min_quality;			/**< Minimum accepted quality. */
//     int max_quality;			/**< Maximum accepted quality. */
//     int base_quality;			/**< Base quality for normalization. */
//     int max_nts_mismatch;		/**< Number of last nts to compute statistics for filtering. */
//     int max_n_per_read;			/**< Number of last nts to compute statistics for filtering. */
//     int min_read_length;		/**< Minimum read length allowed. */
//     int max_read_length;		/**< Maximum read length allowed. */
//     int rtrim_nts;			/**< Number of last nts to compute statistics for preprocessing. */
//     int ltrim_nts;			/**< Number of first nts to compute statistics for preprocessing. */
//     int rfilter_nts;			/**< Number of last nts to compute statistics for filtering. */
//     int lfilter_nts;			/**< Number of first nts to compute statistics for filtering. */
//     int prepro_step;			/**< Flag for preprocessing. */
//     int filter_step;			/**< Flag for filtering. */
//     int qc_step;			/**< Flag for QC calculation. */
//     int num_sources;			/**< Number of sources. */
//     int kmers_on;			/**< Flag for kmers calculation. */
//     int cpu_num_threads;		/**< Number of CPU threads. */
//     int cg_flag;			/**< Chaos game flag. */
//     int k_cg;				/**< Word size k in chaos game for genomic signature. */
//     char* genomic_signature_input;	/**< Genomic signature filename (reference genome). */
//     source_t* source_p;			/**< Fastq sources. */
//     char* report_directory;		/**< Output directory to write HTML reports and data files. */
// } results_server_input_t;

/**
* @brief Writer server parameters
* 
* Structure containing parameters to pass to the writer server
*/
// typedef struct writer_server_input {
//     int rtrim_nts;			/**< Number of last nts to trim in case of preprocessing. */
//     int ltrim_nts;			/**< Number of first nts to trim in case of preprocessing. */
//     int num_sources;			/**< Number of sources. */
//     source_t* source_p;			/**< Fastq sources. */
//     char* output_directory;		/**< Output directory to write the fastq result files. */
// } writer_server_input_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Performs qc, preprocessing and/or filtering of a fastq single end file
*  @param batch_size batch size of fastq reads to load from disk
*  @param max_fastq_batch_list_length maximum length of the batch lists
*  @param num_blocks number of GPU blocks launched
*  @param num_threads number of GPU threads per block launched
*  @param cpu_num_threads number of CPU threads launched for OpenMP paralellized code
*  @param filename fastq filename
*  @param output_directory output directory where output report and files will be written 
*  @param min_quality min quality value accepted
*  @param max_quality max quality value accepted
*  @param base_quality base quality to normalize quality values
*  @param begin_quality_nt initial nucleotide position from which quality will be computed for mean and median quality value calculations
*  @param end_quality_nt final nucleotide position to which quality will be computed for mean and median quality value calculations
*  @param max_nts_mismatch maximum number of nucleotide mismatches for filter purposes
*  @param max_n_per_read maximum number of N positions in a read for filter purposes
*  @param min_read_length mininum read length for filter purposes
*  @param max_read_length maximumum read length for filter purposes 
*  @param rtrim_nts right nucleotides to process for preprocessing purposes 
*  @param ltrim_nts left nucleotides to process for preprocessing purposes 
*  @param rfilter_nts right nucleotides to process for filter purposes 
*  @param lfilter_nts left nucleotides to process for filter purposes
*  @param prepro_step preprocessing flag
*  @param filter_step filtering flag
*  @param qc_step qc flag
*  @param kmers_on kmers flag
*  @param cg_flag chaos game flag
*  @param k_cg word size for the chaos game
*  @param genomic_signature_input filename containing genomic signature of the reference genome
*  @return void
*  
*  Performs qc, preprocessing and/or filtering of a fastq single end file.
*  Inside the functions the needed threads are set up and launched 
*/
void kernel_prepro_fastq_single_end(size_t batch_size, int max_fastq_batch_list_length, int num_blocks, int num_threads, int cpu_num_threads, int cpu_qc_calc_num_threads, char* filename, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input);

/**
*  @brief Performs qc, preprocessing and/or filtering of a fastq paired end file
*  @param batch_size batch size of fastq reads to load from disk
*  @param max_fastq_batch_list_length maximum length of the batch lists
*  @param num_blocks number of GPU blocks launched
*  @param num_threads number of GPU threads per block launched
*  @param cpu_num_threads number of CPU threads launched for OpenMP paralellized code
*  @param filename1 fastq filename paired end 1
*  @param filename2 fastq filename paired end 2 
*  @param output_directory output directory where output report and files will be written 
*  @param min_quality min quality value accepted
*  @param max_quality max quality value accepted
*  @param base_quality base quality to normalize quality values
*  @param begin_quality_nt initial nucleotide position from which quality will be computed for mean and median quality value calculations
*  @param end_quality_nt final nucleotide position to which quality will be computed for mean and median quality value calculations
*  @param max_nts_mismatch maximum number of nucleotide mismatches for filter purposes
*  @param max_n_per_read maximum number of N positions in a read for filter purposes
*  @param min_read_length mininum read length for filter purposes
*  @param max_read_length maximumum read length for filter purposes 
*  @param rtrim_nts right nucleotides to process for preprocessing purposes 
*  @param ltrim_nts left nucleotides to process for preprocessing purposes 
*  @param rfilter_nts right nucleotides to process for filter purposes 
*  @param lfilter_nts left nucleotides to process for filter purposes
*  @param prepro_step preprocessing flag
*  @param filter_step filtering flag
*  @param qc_step qc flag
*  @param kmers_on kmers flag
*  @param cg_flag chaos game flag
*  @param k_cg word size for the chaos game
*  @param genomic_signature_input filename containing genomic signature of the reference genome
*  @return void
*  
*  Performs qc, preprocessing and/or filtering of a fastq paired end file.
*  Inside the functions the needed threads are set up and launched
*/
void kernel_prepro_fastq_paired_end(size_t batch_size, int max_fastq_batch_list_length, int num_blocks, int num_threads, int cpu_num_threads, int cpu_qc_calc_num_threads, char* filename1, char* filename2, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input);

#endif
