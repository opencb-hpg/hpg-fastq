#ifndef OPTIONS_H
#define OPTIONS_H

/*
 * hpg-fastq.h
 *
 *  Created on: Aug 29, 2012
 *      Author: imedina
 */

#include <stdlib.h>
#include <string.h>

#include <argtable2.h>
#include <libconfig.h>

#include <log.h>

#define DEFAULT_MIN_READ_LENGTH		50
#define DEFAULT_MAX_READ_LENGTH		200
#define DEFAULT_MIN_READ_QUALITY		20
#define DEFAULT_MAX_READ_QUALITY		60
#define DEFAULT_MAX_NTS_OUT_QUALITY	3
#define DEFAULT_MAX_N_PER_READ			0
#define DEFAULT_START_QUALITY_NT		0
#define DEFAULT_END_QUALITY_NT			1024

#define DEFAULT_GPU_NUM_BLOCKS			16
#define DEFAULT_GPU_NUM_THREADS		512
#define DEFAULT_GPU_NUM_DEVICES		0
#define DEFAULT_CPU_NUM_THREADS		4
#define DEFAULT_QC_CALC_NUM_THREADS	0
#define DEFAULT_BATCH_SIZE_MB			64*1000000
#define DEFAULT_BATCH_LIST_SIZE		4

#define DEFAULT_K_IN_CHAOS_GAME		10

#define DEFAULT_LOG_LEVEL				2

#define NUM_OPTIONS					30

typedef struct options {
	char *fastq_file; /**< FASTQ file used as input. */
	char *fastq1_file; /**< FASTQ1 file used as input. */
	char *fastq2_file; /**< FASTQ2 file used as input. */
	char *config_file; /**< Path to the configuration file */
	char *genomic_signature_input; /**< Filename template for the main output file. */
	char *output_directory; /**< Directory where the output files will be stored. */

	int min_read_length;	// 50
	int max_read_length;	// 200
	int min_quality;  	// 20
	int max_quality;  	// 60
	int max_nts_out_quality;
	int max_n_per_read;
	int start_quality_nt;
	int end_quality_nt;
	int phred_quality;
	int rtrim_nts;
	int ltrim_nts;
	int rfilter_nts;
	int lfilter_nts;
	int kmers_flag;

	int cg_flag; //Chaos Game flag
	int k_cg;    //7

	// variables to store the value options
	//    struct arg_int *gpu_num_blocks =  DEFAULT_GPU_NUM_BLOCKS;   // 16
	//    struct arg_int *gpu_num_threads =  DEFAULT_GPU_NUM_THREADS;  // 512
	//    struct arg_int *gpu_num_devices =  DEFAULT_GPU_NUM_DEVICES;  // -1
	int cpu_num_threads;  			// 2
	int cpu_qc_calc_num_threads; 	//0
	int batch_size; 				// 64MB
	int batch_list_size;  			// 4

	int log_level; /**< desc */
	char *log_file; /**< desc */
	int verbose; /**< desc */
	int help; /**< desc */
	int time; /**< desc */

} options_t;


options_t *options_new(void);


void options_free(options_t *options);



/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 *
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
void** argtable_options_new(void);


/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 *
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void argtable_options_free(void **argtable_options);


/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 *
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
options_t *read_CLI_options(void **argtable_options, options_t *options);



/**
 * @brief Reads the configuration parameters of the application.
 * @param filename file the options data are read from
 * @param options_data options values (vcf filename, output directory...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 *
 * Reads the basic configuration parameters of the application. If the configuration file can't be
 * read, these parameters should be provided via the command-line interface.
 */
int read_config_file(const char *filename, options_t *options);



options_t *parse_options(int argc, char **argv);


void usage(void **argtable);

#endif

//typedef struct argtable_options {
//	/*	IO options	*/
//	struct arg_file *fastq_file; /**< VCF file used as input. */
//	struct arg_file *fastq1_file; /**< PED file used as input. */
//	struct arg_file *fastq2_file; /**< PED file used as input. */
//	struct arg_file *config_file; /**< Path to the configuration file */
//	struct arg_file *genomic_signature_input; /**< Filename template for the main output file. */
//	struct arg_str *output_directory; /**< Directory where the output files will be stored. */
//
//	//    struct arg_int qc_flag;
//	//    struct arg_int *filter_flag;
//	//    struct arg_int *prepro_flag;
//	struct arg_int *min_read_length;	// 50
//	struct arg_int *max_read_length;	// 200
//	struct arg_int *min_quality;  	// 20
//	struct arg_int *max_quality;  	// 60
//	struct arg_int *max_nts_out_quality;
//	struct arg_int *max_n_per_read;
//	struct arg_int *start_quality_nt;
//	struct arg_int *end_quality_nt;
//	struct arg_int *phred_quality;
//	struct arg_int *rtrim_nts;
//	struct arg_int *ltrim_nts;
//	struct arg_int *rfilter_nts;
//	struct arg_int *lfilter_nts;
//	struct arg_lit *kmers_flag;
//
//	struct arg_lit *cg_flag; //Chaos Game flag
//	struct arg_int *k_cg;    //7
//
//	// variables to store the value options
//	//    struct arg_int *gpu_num_blocks =  DEFAULT_GPU_NUM_BLOCKS;   // 16
//	//    struct arg_int *gpu_num_threads =  DEFAULT_GPU_NUM_THREADS;  // 512
//	//    struct arg_int *gpu_num_devices =  DEFAULT_GPU_NUM_DEVICES;  // -1
//	struct arg_int *cpu_num_threads;  			// 2
//	//    struct arg_int *cpu_qc_calc_num_threads; 	//0
//	struct arg_int *batch_size; 				// 64MB
//	struct arg_int *batch_list_size;  			// 4
//
//	struct arg_int *log_level; /**< desc */
//	struct arg_file *log_file; /**< desc */
//	struct arg_lit *verbose; /**< desc */
//	struct arg_lit *help; /**< desc */
//	struct arg_lit *time; /**< desc */
//
//	int num_options;
//} argtable_options_t;
