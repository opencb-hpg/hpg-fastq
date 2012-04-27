
#ifndef PREPRO_COMMONS_H
#define PREPRO_COMMONS_H

#include "fastq_batch_reader.h"
#include "fastq_read.h"

// #include "chaos_game.h"
// #include "qc_batch.h"
#include "qc_report.h"
// #include "prepro.h"
#include "prepro_batch.h"

/**
* @brief QC calculator server parameters
* 
* Structure containing parameters to pass to the QC calculator server
*/
typedef struct qc_calc_server_input {
    int num_gpu_devices;		/**< Number of GPU devices. */
    int cpu_num_threads;		/**< Number of CPU threads. */
    int gpu_device_id[256];		/**< GPU device identifiers to use. */
    int num_blocks;			/**< Number of GPU blocks to launch. */
    int num_threads;			/**< Number of GPU threads by warp. */
    int min_quality;			/**< Minimum accepted quality. */
    int max_quality;			/**< Maximum accepted quality. */
    int begin_quality_nt;		/**< First nt in the read to compute mean and median quality. */
    int end_quality_nt;			/**< Last nt in the read to compute mean and median quality. */
    int rtrim_nts;			/**< Number of last nts to compute statistics for preprocessing. */
    int ltrim_nts;			/**< Number of first nts to compute statistics for preprocessing. */
    int rfilter_nts;			/**< Number of last nts to compute statistics for filtering. */
    int lfilter_nts;			/**< Number of first nts to compute statistics for filtering. */
    int kmers_on;			/**< Flag for kmers calculation. */
    fastq_batch_list_t* batch_list_p;	/**< Pointer to the fastq batch list. */
} qc_calc_server_input_t;

/**
* @brief Results server parameters
* 
* Structure containing parameters to pass to the results server
*/
typedef struct results_server_input {
    int num_blocks;			/**< Number of GPU blocks launched. */
    int num_threads;			/**< Number of GPU threads launched by warp. */
    int min_quality;			/**< Minimum accepted quality. */
    int max_quality;			/**< Maximum accepted quality. */
    int base_quality;			/**< Base quality for normalization. */
    int max_nts_mismatch;		/**< Number of last nts to compute statistics for filtering. */
    int max_n_per_read;			/**< Number of last nts to compute statistics for filtering. */
    int min_read_length;		/**< Minimum read length allowed. */
    int max_read_length;		/**< Maximum read length allowed. */
    int rtrim_nts;			/**< Number of last nts to compute statistics for preprocessing. */
    int ltrim_nts;			/**< Number of first nts to compute statistics for preprocessing. */
    int rfilter_nts;			/**< Number of last nts to compute statistics for filtering. */
    int lfilter_nts;			/**< Number of first nts to compute statistics for filtering. */
    int prepro_step;			/**< Flag for preprocessing. */
    int filter_step;			/**< Flag for filtering. */
    int qc_step;			/**< Flag for QC calculation. */
    int num_sources;			/**< Number of sources. */
    int kmers_on;			/**< Flag for kmers calculation. */
    int cpu_num_threads;		/**< Number of CPU threads. */
    int cg_flag;			/**< Chaos game flag. */
    int k_cg;				/**< Word size k in chaos game for genomic signature. */
    char* genomic_signature_input;	/**< Genomic signature filename (reference genome). */
    source_t* source_p;			/**< Fastq sources. */
    char* report_directory;		/**< Output directory to write HTML reports and data files. */
} results_server_input_t;

/**
* @brief Writer server parameters
* 
* Structure containing parameters to pass to the writer server
*/
typedef struct writer_server_input {
    int rtrim_nts;			/**< Number of last nts to trim in case of preprocessing. */
    int ltrim_nts;			/**< Number of first nts to trim in case of preprocessing. */
    int num_sources;			/**< Number of sources. */
    source_t* source_p;			/**< Fastq sources. */
    char* output_directory;		/**< Output directory to write the fastq result files. */
} writer_server_input_t;

/* **************************************
 *  		Functions		*
 * *************************************/

void preprocessing_read(int i, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p);
void filtering_read(int i, int length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p);
void qc_per_read(int i, int length, int base_quality, int* mean_read_quality, long* read_sum_length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p, qc_report_t* qc_report, int cpu_num_threads);
void kmers_per_read(int i, int length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p, pthread_mutex_t kmers_lock, int cpu_num_threads);
void kmers_thread_per_read(int i, int length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, int* kmers_valid_thread_counter, int* kmers_invalid_thread_counter, status_batch_t* status_batch_p, int cpu_num_threads);
void qc_per_batch(int source_id, results_server_input_t* input_p, qc_batch_t* qc_batch_p, qc_report_t* qc_report, int valid);
void qc_results(results_server_input_t* input_p, long* read_sum_length, qc_report_t* qc_report, int valid);
void preprocessing_read(int i, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p);
void filtering_read(int i, int read_length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p);
int kmers_index(char* nt);

#endif /* PREPRO_COMMONS_H */
