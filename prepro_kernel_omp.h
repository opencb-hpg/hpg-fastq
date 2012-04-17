
#ifndef PREPRO_KERNEL_OMP_H
#define PREPRO_KERNEL_OMP_H

/**
*  @brief OpenMP Kernel implementation for qc, preprocessing and filtering
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
void cpu_prepro_trim(const int *d_read_batch_data_indices_p, char* d_seq_data_p, char* d_quality_data_p, qc_read_t* d_gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int begin_quality_nt, int end_quality_nt);

/**
*  @brief OpenMP Kernel Kernel for qc, preprocessing and filtering
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
void cpu_prepro_filter_trim(const int *read_batch_data_indices_p, char* d_seq_data_p, char* d_quality_data_p, qc_read_t* gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int begin_quality_nt, int end_quality_nt);

#endif	/*  PREPRO_KERNEL_OMP_H	 */

