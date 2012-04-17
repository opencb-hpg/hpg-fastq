
#include "commons.h"
#include "cuda_commons.h"
#include "fastq_file.h"
#include "fastq_batch_list.h"
#include "fastq_batch_reader.h"
#include "fastq_read.h"
#include "log.h"
#include "prepro_kernel_omp.h"
#include "qc_batch.h"
#include "qc_report.h"

//---------------------------------------------------------------------------------------------
//  CPU kernel for preprocessing/filtering rtrim_nts == rfilter_nts and ltrim_nts == lfilter_nts
//  both sides quality is calculated only once (CPU implementation)
//---------------------------------------------------------------------------------------------

void cpu_prepro_trim(const int *d_read_batch_data_indices_p, char* d_seq_data_p, char* d_quality_data_p, qc_read_t* d_gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int begin_quality_nt, int end_quality_nt) {
	char nt, quality;
	qc_read_t* qc_read_p;
	int i, j, k, read_length, nt_position;
	int rtrim_nt_position, ltrim_nt_position, begin_nt_position, end_nt_position;
	int rtrim_nts_quality, ltrim_nts_quality;
	int trim_read_quality, ltrim_read_quality, rtrim_read_quality;
	int accumulated_total_quality = 0, accumulated_mean_quality = 0;
	int out_of_quality = 0, left_out_of_quality = 0, right_out_of_quality = 0;
	int rtrim_nts_accumulated_quality = 0, ltrim_nts_accumulated_quality = 0;
	int mean_quality, median_quality = 0, median_aux_count = 0;
	int quality_window, half_quality_window;
	int counters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int lcounters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int rcounters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int quality_histogram[MAX_QUALITY_VALUE + PHRED64];
	
	for (k=0; k < num_reads; k++) {
		qc_read_p = &d_gpu_result_p[k];

		i = d_read_batch_data_indices_p[k];
		read_length = d_read_batch_data_indices_p[k+1] - i;

		if (begin_quality_nt >= read_length) { begin_quality_nt = 1;}
		if (end_quality_nt > read_length) { end_quality_nt = read_length;}
		quality_window = end_quality_nt - begin_quality_nt;
		half_quality_window = (quality_window + 1) / 2;

		j = i;
		rtrim_nt_position = j + read_length - rtrim_nts;
		ltrim_nt_position = j + ltrim_nts;
		begin_nt_position = j + begin_quality_nt;
		end_nt_position = j + end_quality_nt;
		memset(&quality_histogram, 0, sizeof(quality_histogram));

		while((nt = d_seq_data_p[i++]) != '\0') {
		  quality = d_quality_data_p[j++];
		  
		  nt_position = (nt & 7);
		  counters[nt_position]++;

		  accumulated_total_quality += quality;

		  if ((j >= begin_nt_position) && (j <=  end_nt_position)) {	// we look only in the nt window of interest
		    accumulated_mean_quality += quality;	// for mean calculation
		    quality_histogram[quality]++;		// for median calculation

		    if ((quality < min_quality) || (quality > max_quality)) {
		      out_of_quality++;			//for out-of-range qualities
		    }
		  }

		  if (j >= rtrim_nt_position) {
		    rtrim_nts_accumulated_quality += quality;

		    if ((quality < min_quality) || (quality > max_quality)) {
		      right_out_of_quality++;
		    }

		    rcounters[nt_position]++;
		  }

		  if (j <= ltrim_nt_position) {
		    ltrim_nts_accumulated_quality += quality;

		    if ((quality < min_quality) || (quality > max_quality)) {
		      left_out_of_quality++;
		    }

		    lcounters[nt_position]++;
		  }
		}

		median_aux_count = 0;

		for (int h=0; h < MAX_QUALITY_VALUE + PHRED64; h++) {
		  median_aux_count += quality_histogram[h];

		  if (median_aux_count >= half_quality_window) {
		    median_quality = h;
		    break;
		  }
		}

		trim_read_quality = accumulated_total_quality - rtrim_nts_accumulated_quality - ltrim_nts_accumulated_quality;
		ltrim_read_quality = accumulated_total_quality - ltrim_nts_accumulated_quality;
		rtrim_read_quality = accumulated_total_quality - rtrim_nts_accumulated_quality;

		qc_read_p->counters[A] = (lcounters[A] << 22) + (rcounters[A] << 12) + counters[A];
		qc_read_p->counters[C] = (lcounters[C] << 22) + (rcounters[C] << 12) + counters[C];
		qc_read_p->counters[G] = (lcounters[G] << 22) + (rcounters[G] << 12) + counters[G];
		qc_read_p->counters[T] = (lcounters[T] << 22) + (rcounters[T] << 12) + counters[T];
		qc_read_p->counters[N] = (lcounters[N] << 22) + (rcounters[N] << 12) + counters[N];
		qc_read_p->counters[TOTAL] = read_length - 1;

		if (accumulated_mean_quality == accumulated_total_quality) {
		  qc_read_p->counters[MEAN_MEDIAN_QUALITY] = (median_quality << 8) + (accumulated_total_quality / (read_length - 1));
		} else {
		  end_quality_nt = max(read_length, end_quality_nt);
		  qc_read_p->counters[MEAN_MEDIAN_QUALITY] = (median_quality << 8) + (accumulated_mean_quality / (quality_window + 1));
		}

		mean_quality = GET_0_8(qc_read_p->counters[MEAN_MEDIAN_QUALITY]);
		rtrim_nts_quality = (rtrim_nts == 0) ? mean_quality : rtrim_nts_accumulated_quality / rtrim_nts;
		ltrim_nts_quality = (ltrim_nts == 0) ? mean_quality : ltrim_nts_accumulated_quality / ltrim_nts;

		qc_read_p->counters[SIDE_NTS_QUALITY] = (ltrim_nts_quality << 24) + (ltrim_nts_quality << 16) + (rtrim_nts_quality << 8) + rtrim_nts_quality;

		trim_read_quality = (trim_read_quality / (read_length - rtrim_nts - ltrim_nts - 1));
		ltrim_read_quality = (ltrim_read_quality / (read_length - ltrim_nts - 1));
		rtrim_read_quality = (rtrim_read_quality / (read_length - rtrim_nts - 1));

		qc_read_p->counters[TRIM_QUALITY] = (ltrim_read_quality << 24) + (rtrim_read_quality << 16) + (trim_read_quality << 8) + (accumulated_total_quality / (read_length - 1));
		qc_read_p->counters[OUT_OF_QUALITY] = (left_out_of_quality << 22) + (right_out_of_quality << 12) + out_of_quality;

		// clean values for next read iteration
		accumulated_total_quality = 0;
		accumulated_mean_quality = 0;
		out_of_quality = 0;
		left_out_of_quality = 0;
		right_out_of_quality = 0;
		rtrim_nts_accumulated_quality = 0;
		ltrim_nts_accumulated_quality = 0;
		mean_quality = 0;
		median_quality = 0;
		median_aux_count = 0;

		for (int a=0; a < COUNTERS_SIZE; a++) { 
		  counters[a] = 0;
		  lcounters[a] = 0;
		  rcounters[a] = 0;
		}
	}
}

//---------------------------------------------------------------------------------------------
//  CPU kernel for preprocessing/filtering rtrim_nts != rfilter_nts and ltrim_nts != lfilter_nts
//  trim and filter quality is calculated separately (CPU implementation)
//---------------------------------------------------------------------------------------------

void cpu_prepro_filter_trim(const int *read_batch_data_indices_p, char* d_seq_data_p, char* d_quality_data_p, qc_read_t* gpu_result_p, int num_reads, int min_quality, int max_quality, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int begin_quality_nt, int end_quality_nt) {
	char nt, quality;
	qc_read_t* qc_read_p;
	int i, j, k, read_length, nt_position;
	int rtrim_nt_position, ltrim_nt_position, rfilter_nt_position, lfilter_nt_position, begin_nt_position, end_nt_position;
	int rtrim_nts_quality, ltrim_nts_quality, rfilter_nts_quality, lfilter_nts_quality;
	int trim_read_quality, ltrim_read_quality, rtrim_read_quality;
	int lfilter_read_quality, rfilter_read_quality;
	int accumulated_total_quality = 0, accumulated_mean_quality = 0;
	int out_of_quality = 0, left_out_of_quality = 0, right_out_of_quality = 0;
	int rtrim_nts_accumulated_quality = 0, ltrim_nts_accumulated_quality = 0;
	int rfilter_nts_accumulated_quality = 0, lfilter_nts_accumulated_quality = 0;
	int mean_quality, median_quality = 0, median_aux_count = 0;
	int quality_window = end_quality_nt - begin_quality_nt;
	int half_quality_window = (quality_window + 1) / 2;
	int counters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int lcounters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int rcounters[COUNTERS_SIZE] = {0,0,0,0,0,0,0,0,0,0};
	int quality_histogram[MAX_QUALITY_VALUE + PHRED64];
	memset(&quality_histogram, 0, sizeof(quality_histogram));

	for (k=0; k < num_reads; k++) {
		qc_read_p = &gpu_result_p[k];

		i = read_batch_data_indices_p[k];
		read_length = read_batch_data_indices_p[k+1] - i;
		j = i;
		rtrim_nt_position = j + read_length - rtrim_nts;
		ltrim_nt_position = j + ltrim_nts;
		rfilter_nt_position = j + read_length - rfilter_nts;
		lfilter_nt_position = j + lfilter_nts;
		begin_nt_position = j + begin_quality_nt;
		end_nt_position = j + end_quality_nt;

		while ((nt = d_seq_data_p[i++]) != '\0') {
		  quality = d_quality_data_p[j++];  
		  
		  nt_position = (nt & 7);
		  counters[nt_position]++;

		  accumulated_total_quality += quality;

		  if ((j >= begin_nt_position) && (j <=  end_nt_position)) {	// we look only in the nt window of interest
		    accumulated_mean_quality += quality;	// for mean calculation
		    quality_histogram[quality]++;		// for median calculation

		    if ((quality < min_quality) || (quality > max_quality)) {
		      out_of_quality++;			//for out-of-range qualities
		    }
		  }

		  if (j >= rtrim_nt_position) {
		    rtrim_nts_accumulated_quality += quality;

		    if ((quality < min_quality) || (quality > max_quality)) {
		      right_out_of_quality++;
		    }

		    rcounters[nt_position]++;
		  }

		  if (j <= ltrim_nt_position) {
		    ltrim_nts_accumulated_quality += quality;

		    if ((quality < min_quality) || (quality > max_quality)) {
		      left_out_of_quality++;
		    }

		    lcounters[nt_position]++;
		  }

		  if (j >= rfilter_nt_position) {
		    rfilter_nts_accumulated_quality += quality;
		  }

		  if (j <= lfilter_nt_position) {
		    lfilter_nts_accumulated_quality += quality;
		  }
		}

		median_aux_count = 0;

		for (int h=0; h < MAX_QUALITY_VALUE + PHRED64; h++) {
		  median_aux_count += quality_histogram[h];

		  if (median_aux_count >= half_quality_window) {
		    median_quality = h;
		    break;
		  }
		}

		trim_read_quality = accumulated_total_quality - rtrim_nts_accumulated_quality - ltrim_nts_accumulated_quality;
		ltrim_read_quality = accumulated_total_quality - ltrim_nts_accumulated_quality;
		rtrim_read_quality = accumulated_total_quality - rtrim_nts_accumulated_quality;
		lfilter_read_quality = accumulated_total_quality - lfilter_nts_accumulated_quality;
		rfilter_read_quality = accumulated_total_quality - rfilter_nts_accumulated_quality;

		qc_read_p->counters[A] = (lcounters[A] << 22) + (rcounters[A] << 12) + counters[A];
		qc_read_p->counters[C] = (lcounters[C] << 22) + (rcounters[C] << 12) + counters[C];
		qc_read_p->counters[G] = (lcounters[G] << 22) + (rcounters[G] << 12) + counters[G];
		qc_read_p->counters[T] = (lcounters[T] << 22) + (rcounters[T] << 12) + counters[T];
		qc_read_p->counters[N] = (lcounters[N] << 22) + (rcounters[N] << 12) + counters[N];
		qc_read_p->counters[TOTAL] = read_length - 1;

		if (accumulated_mean_quality == accumulated_total_quality) {
		  qc_read_p->counters[MEAN_MEDIAN_QUALITY] = (median_quality << 8) + (accumulated_total_quality / (read_length - 1));
		} else {
		  end_quality_nt = max(read_length, end_quality_nt);
		  qc_read_p->counters[MEAN_MEDIAN_QUALITY] = (median_quality << 8) + (accumulated_mean_quality / (quality_window + 1));
		}

		mean_quality = GET_0_8(qc_read_p->counters[MEAN_MEDIAN_QUALITY]);
		rtrim_nts_quality = (rtrim_nts == 0) ? mean_quality : rtrim_nts_accumulated_quality / rtrim_nts;
		ltrim_nts_quality = (ltrim_nts == 0) ? mean_quality : ltrim_nts_accumulated_quality / ltrim_nts;
		rfilter_nts_quality = (rfilter_nts == 0) ? mean_quality : rfilter_nts_accumulated_quality / rfilter_nts;
		lfilter_nts_quality = (lfilter_nts == 0) ? mean_quality : lfilter_nts_accumulated_quality / lfilter_nts;

		qc_read_p->counters[SIDE_NTS_QUALITY] = (ltrim_nts_quality << 24) + (lfilter_nts_quality << 16) + (rfilter_nts_quality << 8) + rtrim_nts_quality;

		trim_read_quality = (trim_read_quality / (read_length - rtrim_nts - ltrim_nts - 1));
		ltrim_read_quality = (ltrim_read_quality / (read_length - ltrim_nts - 1));
		rtrim_read_quality = (rtrim_read_quality / (read_length - rtrim_nts - 1));

		qc_read_p->counters[TRIM_QUALITY] = (ltrim_read_quality << 24) + (rtrim_read_quality << 16) + (trim_read_quality << 8) + (accumulated_total_quality / (read_length - 1));
		qc_read_p->counters[OUT_OF_QUALITY] = (left_out_of_quality << 22) + (right_out_of_quality << 12) + out_of_quality;

		// clean values for next read iteration
		accumulated_total_quality = 0;
		accumulated_mean_quality = 0;
		out_of_quality = 0;
		left_out_of_quality = 0;
		right_out_of_quality = 0;
		rtrim_nts_accumulated_quality = 0;
		ltrim_nts_accumulated_quality = 0;
		mean_quality = 0;
		median_quality = 0;
		median_aux_count = 0;

		for (int a=0; a < COUNTERS_SIZE; a++) { 
		  counters[a] = 0;
		  lcounters[a] = 0;
		  rcounters[a] = 0;
		}
	}
}





