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


#ifndef PREPRO_KERNEL_CUDA_H
#define PREPRO_KERNEL_CUDA_H

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief CUDA Kernel for qc, preprocessing and filtering
*  @param d_read_batch_data_p pointer to data 
*  @param d_read_batch_data_indices_p pointer to data_indices
*  @param[out] d_gpu_result_p pointer to results structure
*  @param num_reads number of reads in data vector
*  @param min_quality min accepted quality
*  @param max_quality max accepted quality
*  @param rtrim_nts number of last nts to compute statistics for preprocessing and filtering (near 3')
*  @param ltrim_nts number of first nts to compute statistics for preprocessing and filtering (near 5')
*  @param begin_quality_nt first nt in the read to compute mean and median quality
*  @param end_quality_nt last nt in the read to compute mean and median quality
*  @return void 
*  
*  Performs qc, preprocessing and/or filtering calculations. All scenaries calculations are returned
*  and the results thread decides which reads are filtered or preprocessed.
*  Number of nucleotides of filtering and preprocessing are the SAME in this kernel implementation.
*/
__global__ void kernel_prepro(int *d_read_batch_data_indices_p, char* d_seq_p, char* d_quality_p, qc_read_t* d_gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int begin_quality_nt, int end_quality_nt);

/**
*  @brief CUDA Kernel for qc, preprocessing and filtering
*  @param d_read_batch_data_p pointer to data 
*  @param d_read_batch_data_indices_p pointer to data_indices
*  @param[out] d_gpu_result_p pointer to results structure
*  @param num_reads number of reads in data vector
*  @param min_quality min accepted quality
*  @param max_quality max accepted quality
*  @param rtrim_nts number of last nts to compute statistics for preprocessing (near 3')
*  @param ltrim_nts number of first nts to compute statistics for preprocessing (near 5')
*  @param rfilter_nts number of last nts to compute statistics for filtering (near 3')
*  @param lfilter_nts number of first nts to compute statistics for filtering (near 5')
*  @param begin_quality_nt first nt in the read to compute mean and median quality
*  @param end_quality_nt last nt in the read to compute mean and median quality
*  @return void 
*  
*  Performs qc, preprocessing and/or filtering calculations. All scenaries calculations are returned
*  and the results thread decides which reads are filtered or preprocessed.
*  Number of nucleotides of filtering and preprocessing are DIFFERENT in this kernel implementation.
*/
__global__ void kernel_prepro(int *d_read_batch_data_indices_p, char* d_seq_p, char* d_quality_p, qc_read_t* d_gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int begin_quality_nt, int end_quality_nt);

/**
*  @brief CUDA Kernel for grep operation
*  @param d_read_batch_data_p pointer to data 
*  @param d_read_batch_data_indices_p pointer to data_indices
*  @param d_match_p matching string
*  @param match_length length of the matching string
*  @param num_reads number of reads in the data vector
*  @return void 
*  
*  Performs grep operation. Finds a string match in the data vector.
*/
__global__ void kernel_grep(char *d_read_batch_data_p, int *d_read_batch_data_indices_p, char *d_match_p, int match_length, int num_reads);

/**
*  @brief CUDA Kernel for qc
*  @param d_read_batch_data_indices_p pointer to sequence and quality indices
*  @param d_seq_p pointer to sequences vector
*  @param d_quality_p pointer to qualities vector
*  @param[out] d_gpu_result_p pointer to results structure
*  @param num_reads number of reads in data vector
*  @param begin_quality_nt first nt in the read to compute mean and median quality
*  @param end_quality_nt last nt in the read to compute mean and median quality
*  @return void 
*  
*  Performs qc calculations.
*/
__global__ void kernel_qc(int *d_read_batch_data_indices_p, char* d_seq_p, char* d_quality_p, qc_read_t* d_gpu_result_p, int num_reads, int begin_quality_nt, int end_quality_nt);

#endif	/*	PREPRO_KERNEL_CUDA_H	*/
