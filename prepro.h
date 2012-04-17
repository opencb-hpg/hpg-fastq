
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
#include "cuda_commons.h"
#include "list.h"
#include "log.h"
#include "file_utils.h"
#include "system_utils.h"

#define VALID_READS_FILE_SUFFIX 	".valid"
#define INVALID_READS_FILE_SUFFIX 	".invalid"

#define	ALL_READS			1
#define	ONLY_VALID_READS		2
#define	ONLY_INVALID_READS		3

/**
*  @brief Performs qc, preprocessing and/or filtering of a fastq single end file
*  @param batch_size batch size of fastq reads to load from disk
*  @param max_fastq_batch_list_length maximum length of the batch lists
*  @param nb_blocks number of GPU blocks launched
*  @param nb_threads number of GPU threads per block launched
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
void kernel_prepro_fastq_single_end(size_t batch_size, int max_fastq_batch_list_length, int nb_blocks, int nb_threads, int cpu_num_threads, char* filename, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input);

/**
*  @brief Performs qc, preprocessing and/or filtering of a fastq paired end file
*  @param batch_size batch size of fastq reads to load from disk
*  @param max_fastq_batch_list_length maximum length of the batch lists
*  @param nb_blocks number of GPU blocks launched
*  @param nb_threads number of GPU threads per block launched
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
void kernel_prepro_fastq_paired_end(size_t batch_size, int max_fastq_batch_list_length, int nb_blocks, int nb_threads, int cpu_num_threads, char* filename1, char* filename2, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input);

#endif