#include "options.h"

#include <argtable2.h>


options_t *options_new(void) {
	options_t *options = (options_t*) calloc (1, sizeof(options_t));

	options->min_read_length = DEFAULT_MIN_READ_LENGTH;
	options->max_read_length = DEFAULT_MAX_READ_LENGTH;
	options->min_quality = DEFAULT_MIN_READ_QUALITY;
	options->max_quality = DEFAULT_MAX_READ_QUALITY;
	options->max_nts_out_quality = DEFAULT_MAX_NTS_OUT_QUALITY;
	options->max_n_per_read = DEFAULT_MAX_N_PER_READ;
	options->start_quality_nt = DEFAULT_START_QUALITY_NT;
	options->end_quality_nt = DEFAULT_END_QUALITY_NT;
	options->kmers_flag = 0;

	options->cg_flag = 0;
	options->k_cg = DEFAULT_K_IN_CHAOS_GAME;

	options->cpu_num_threads = DEFAULT_CPU_NUM_THREADS;
	options->cpu_qc_calc_num_threads = 2;
	options->batch_size = DEFAULT_BATCH_SIZE_MB;
	options->batch_list_size = DEFAULT_BATCH_LIST_SIZE;

	options->log_level = DEFAULT_LOG_LEVEL;
	options->verbose = 0;
	options->time = 0;

	return options;
}


void options_free(options_t *options) {
	if(options == NULL) {
		return;
	}
	if (options->fastq_file)				{ free(options->fastq_file); }
	if (options->fastq1_file)				{ free(options->fastq1_file); }
	if (options->fastq2_file)				{ free(options->fastq2_file); }
	if (options->config_file)				{ free(options->config_file); }
	if (options->genomic_signature_input)	{ free(options->genomic_signature_input); }
	if (options->output_directory)			{ free(options->output_directory); }
	if (options->log_file)					{ free(options->log_file); }
	free(options);
}


void** argtable_options_new(void) {
	void **argtable = (void**)malloc((NUM_OPTIONS+1) * sizeof(void*));	// NUM_OPTIONS +1 to allocate end structure

	// NOTICE that order cannot be changed as is accessed by index in other functions
	argtable[0] = arg_file0("f", "fq,fastq", NULL, "PED file used as input");
	argtable[1] = arg_file0("1", "fq1,fastq1", NULL, "PED file used as input");
	argtable[2] = arg_file0("2", "fq2,fastq2", NULL, "PED file used as input");
	argtable[3] = arg_file0("c", "conf", NULL, "PED file used as input");
	argtable[4] = arg_file0(NULL, "gs-file", NULL, "PED file used as input");
	argtable[5] = arg_str0("o", "outdir", NULL, "PED file used as input");

	argtable[6] = arg_int0("l", "min-read-length", NULL, "Maximum number of batches stored at the same time");
	argtable[7] = arg_int0("L", "max-read-length", NULL, "Maximum number of batches stored at the same time");
	argtable[8] = arg_int0("q", "min-quality", NULL, "Maximum number of batches stored at the same time");
	argtable[9] = arg_int0("Q", "max-quality", NULL, "Maximum number of batches stored at the same time");
	argtable[10] = arg_int0(NULL, "max-nts-out-quality", NULL, "Maximum number of batches stored at the same time");
	argtable[11] = arg_int0("n", "max-n-per-read,max-N-per-read", NULL, "Maximum number of batches stored at the same time");
	argtable[12] = arg_int0(NULL, "start-quality-nt", NULL, "Maximum number of batches stored at the same time");
	argtable[13] = arg_int0(NULL, "end-quality-nt", NULL, "Maximum number of batches stored at the same time");
	argtable[14] = arg_int0("p", "phred-quality", NULL, "Maximum number of batches stored at the same time");
	argtable[15] = arg_int0(NULL, "rtrim-nts", NULL, "Maximum number of batches stored at the same time");
	argtable[16] = arg_int0(NULL, "ltrim-nts", NULL, "Maximum number of batches stored at the same time");
	argtable[17] = arg_int0(NULL, "rfilter-nts", NULL, "Maximum number of batches stored at the same time");
	argtable[18] = arg_int0(NULL, "lfilter-nts", NULL, "Maximum number of batches stored at the same time");
	argtable[19] = arg_lit0("k", "kmers", "Maximum number of batches stored at the same time");

	argtable[20] = arg_lit0(NULL, "cg", "Maximum number of batches stored at the same time");
	argtable[21] = arg_int0(NULL, "k-cg", NULL, "Maximum number of batches stored at the same time");

	argtable[22] = arg_int0("T", "num-threads", NULL, "Maximum number of batches stored at the same time");
	argtable[23] = arg_int0(NULL, "batch-size", NULL, "Maximum number of batches stored at the same time");
	argtable[24] = arg_int0(NULL, "batch-list-size", NULL, "Maximum number of batches stored at the same time");

	argtable[25] = arg_int0(NULL, "log-level", NULL, "Maximum number of batches stored at the same time");
	argtable[26] = arg_file0(NULL, "log-file", NULL, "PED file used as input");
	argtable[27] = arg_lit0("v", "verbose", "Maximum number of batches stored at the same time");
	argtable[28] = arg_lit0("h", "help", "Maximum number of batches stored at the same time");
	argtable[29] = arg_lit0("t", "time", "Maximum number of batches stored at the same time");

	argtable[30] = arg_end(20);

	return argtable;
}


void argtable_options_free(void **argtable) {
	if(argtable != NULL) {
		arg_freetable(argtable, NUM_OPTIONS+1);	// struct end must also be freed
		free(argtable);
	}
}



int read_config_file(const char *filename, options_t *options) {
	if (filename == NULL || options == NULL) {
		return -1;
	}

	config_t *config = (config_t*) calloc (1, sizeof(config_t));
	int ret_code = config_read_file(config, filename);
	if (ret_code == CONFIG_FALSE) {
		LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
		return -1;
	}

	const char *tmp_string;
	long tmp_int;

	if(config_lookup_string(config, "app.outdir", &tmp_string)) { options->output_directory = strdup(tmp_string); }
	if(config_lookup_int(config, "app.cpu-num-threads", &tmp_int)) { options->cpu_num_threads = (int)tmp_int; }


	config_destroy(config);
	free(config);
//	free(tmp_string);

	return ret_code;
}


/**
 * @brief Initializes an options_t structure from argtable parsed CLI with default values. Notice that options are order dependent.
 * @return A new options_t structure initialized with default values.
 *
 * Initializes the only default options from options_t.
 */
options_t *read_CLI_options(void **argtable, options_t *options) {
	//	options_t *options = (options_t*) calloc (1, sizeof(options_t));

	if (((struct arg_file*)argtable[0])->count) { options->fastq_file = strdup(*(((struct arg_file*)argtable[0])->filename)); }
	if (((struct arg_file*)argtable[1])->count) { options->fastq1_file = strdup(*(((struct arg_file*)argtable[1])->filename)); }
	if (((struct arg_file*)argtable[2])->count) { options->fastq2_file = strdup(*(((struct arg_file*)argtable[2])->filename)); }
	if (((struct arg_file*)argtable[3])->count) { options->config_file = strdup(*(((struct arg_file*)argtable[3])->filename)); }
	if (((struct arg_file*)argtable[4])->count) { options->genomic_signature_input = strdup(*(((struct arg_file*)argtable[4])->filename)); }
	if (((struct arg_str*)argtable[5])->count) { free(options->output_directory); options->output_directory = strdup(*(((struct arg_str*)argtable[5])->sval)); }

	if (((struct arg_int*)argtable[6])->count) { options->min_read_length = *(((struct arg_int*)argtable[6])->ival); }
	if (((struct arg_int*)argtable[7])->count) { options->max_read_length = *(((struct arg_int*)argtable[7])->ival); }
	if (((struct arg_int*)argtable[8])->count) { options->min_quality = *(((struct arg_int *)argtable[8])->ival); }
	if (((struct arg_int*)argtable[9])->count) { options->max_quality = *(((struct arg_int*)argtable[9])->ival); }
	if (((struct arg_int*)argtable[10])->count) { options->max_nts_out_quality = *(((struct arg_int*)argtable[10])->ival); }
	if (((struct arg_int*)argtable[11])->count) { options->max_n_per_read = *(((struct arg_int*)argtable[11])->ival); }
	if (((struct arg_int*)argtable[12])->count) { options->start_quality_nt = *(((struct arg_int*)argtable[12])->ival); }
	if (((struct arg_int*)argtable[13])->count) { options->end_quality_nt = *(((struct arg_int*)argtable[13])->ival); }
	if (((struct arg_int*)argtable[14])->count) { options->phred_quality = *(((struct arg_int*)argtable[14])->ival); }
	if (((struct arg_int*)argtable[15])->count) { options->rtrim_nts = *(((struct arg_int*)argtable[15])->ival); }
	if (((struct arg_int*)argtable[16])->count) { options->ltrim_nts = *(((struct arg_int*)argtable[16])->ival); }
	if (((struct arg_int*)argtable[17])->count) { options->rfilter_nts = *(((struct arg_int*)argtable[17])->ival); }
	if (((struct arg_int*)argtable[18])->count) { options->lfilter_nts = *(((struct arg_int*)argtable[18])->ival); }
	if (((struct arg_int*)argtable[19])->count) { options->kmers_flag = ((struct arg_int*)argtable[19])->count; }

	if (((struct arg_int*)argtable[20])->count) { options->cg_flag = ((struct arg_int*)argtable[20])->count; }
	if (((struct arg_int*)argtable[21])->count) { options->k_cg = *(((struct arg_int*)argtable[21])->ival); }

	if (((struct arg_int*)argtable[22])->count) {options->cpu_num_threads = *(((struct arg_int*)argtable[22])->ival); }
	if (((struct arg_int*)argtable[23])->count) {options->batch_size = *(((struct arg_int*)argtable[23])->ival); }
	if (((struct arg_int*)argtable[24])->count) {options->batch_list_size = *(((struct arg_int*)argtable[24])->ival); }

	if (((struct arg_int*)argtable[25])->count) { options->log_level = *(((struct arg_int*)argtable[25])->ival); }
	if (((struct arg_file*)argtable[26])->count) { options->log_file = strdup(*(((struct arg_file*)argtable[26])->filename)); }
	if (((struct arg_int*)argtable[27])->count) { options->verbose = ((struct arg_int*)argtable[27])->count; }
	if (((struct arg_int*)argtable[28])->count) { options->help = ((struct arg_int*)argtable[28])->count; }
	if (((struct arg_int*)argtable[29])->count) { options->time = ((struct arg_int*)argtable[29])->count; }

	return options;
}


options_t *parse_options(int argc, char **argv) {
	void **argtable = argtable_options_new();
	//	struct arg_end *end = arg_end(10);
	//	void **argtable = argtable_options_get(argtable_options, end);

	options_t *options = options_new();
	if (argc < 2) {
		usage(argtable);
		//exit(0);
	}else {
		int num_errors = arg_parse(argc, argv, argtable);

		if (num_errors > 0) {
			arg_print_errors(stdout, argtable[NUM_OPTIONS], "hpg-fastq");	// struct end is always allocated in the last position
			usage(argtable);
		}else {
			// Check if 'help' option has been provided.
			if(((struct arg_file*)argtable[28])->count > 0) {
				usage(argtable);
				argtable_options_free(argtable);
				options_free(options);
				exit(0);
			}
			// Check if a config file has been provided.
			// CLI options have more priority than config file so this must be read before.
			// --conf option is allocated at index 3
			if(((struct arg_file*)argtable[3])->count > 0) {
				read_config_file(*(((struct arg_file*)argtable[3])->filename), options);
			}

			// Always must be done after arg_parse is called and read_config_file
			options = read_CLI_options(argtable, options);
		}

	}

//	exit:
	argtable_options_free(argtable);
	//	free(end);
	//	free(argtable_options);
	return options;
}

void usage(void **argtable) {
	printf("Usage:\nhpg-fastq {qc | filter | prepro}");
	arg_print_syntaxv(stdout, argtable, "\n");
	arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
}


//argtable_options_t *argtable_options_new(void) {
//	argtable_options_t *argtable_options = (argtable_options_t*) calloc (1, sizeof(argtable_options_t));
//
//	argtable_options->fastq_file = arg_file0("f", "fq,fastq", NULL, "PED file used as input");
//	argtable_options->fastq1_file = arg_file0("1", "fq1,fastq1", NULL, "PED file used as input");
//	argtable_options->fastq2_file = arg_file0("2", "fq2,fastq2", NULL, "PED file used as input");
//	argtable_options->config_file = arg_file0("c", "conf", NULL, "PED file used as input");
//	argtable_options->genomic_signature_input = arg_file0(NULL, "gs-file", NULL, "PED file used as input");
//	argtable_options->output_directory = arg_str1("o", "outdir", NULL, "PED file used as input");
//
//	argtable_options->min_read_length = arg_int0("l", "min-read-length", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->max_read_length = arg_int0("L", "max-read-length", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->min_quality = arg_int0("q", "min-quality", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->max_quality = arg_int0("Q", "max-quality", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->max_nts_out_quality = arg_int0(NULL, "max-nts-out-quality", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->max_n_per_read = arg_int0("n", "max-n-per-read,max-N-per-read", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->start_quality_nt = arg_int0(NULL, "start-quality-nt", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->end_quality_nt = arg_int0(NULL, "end-quality-nt", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->phred_quality = arg_int0("p", "phred-quality", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->rtrim_nts = arg_int0(NULL, "rtrim-nts", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->ltrim_nts = arg_int0(NULL, "ltrim-nts", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->rfilter_nts = arg_int0(NULL, "rfilter-nts", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->lfilter_nts = arg_int0(NULL, "lfilter-nts", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->kmers_flag = arg_lit0("k", "kmers", "Maximum number of batches stored at the same time");
//
//	argtable_options->cg_flag = arg_lit0(NULL, "cg", "Maximum number of batches stored at the same time");
//	argtable_options->k_cg = arg_int0(NULL, "k-cg", NULL, "Maximum number of batches stored at the same time");
//
//	argtable_options->cpu_num_threads = arg_int0("T", "num-threads", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->batch_size = arg_int0(NULL, "batch-size", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->batch_list_size = arg_int0(NULL, "batch-list-size", NULL, "Maximum number of batches stored at the same time");
//
//	argtable_options->log_level = arg_int0(NULL, "log-level", NULL, "Maximum number of batches stored at the same time");
//	argtable_options->log_file = arg_file0(NULL, "log-file", NULL, "PED file used as input");
//	argtable_options->verbose = arg_lit0("v", "verbose", "Maximum number of batches stored at the same time");
//	argtable_options->help = arg_lit0("h", "help", "Maximum number of batches stored at the same time");
//	argtable_options->time = arg_lit0("t", "verbose", "Maximum number of batches stored at the same time");
//
//	argtable_options->num_options = 30;
//
//	return argtable_options;
//}

//void **argtable_options_get(argtable_options_t *argtable_options, struct arg_end *end) {
//	//	struct arg_end *end = arg_end(argtable_options->num_options);
//	void **argtable = (void**)malloc((argtable_options->num_options+1) * sizeof(void*));
//	argtable[0] = argtable_options->fastq_file;
//	argtable[1] = argtable_options->fastq1_file;
//	argtable[2] = argtable_options->fastq2_file;
//	argtable[3] = argtable_options->config_file;
//	argtable[4] = argtable_options->genomic_signature_input;
//	argtable[5] = argtable_options->output_directory;
//
//	argtable[6] = argtable_options->min_read_length;
//	argtable[7] = argtable_options->max_read_length;
//	argtable[8] = argtable_options->min_quality;
//	argtable[9] = argtable_options->max_quality;
//	argtable[10] = argtable_options->max_nts_out_quality;
//	argtable[11] = argtable_options->max_n_per_read;
//	argtable[12] = argtable_options->start_quality_nt;
//	argtable[13] = argtable_options->end_quality_nt;
//	argtable[14] = argtable_options->phred_quality;
//	argtable[15] = argtable_options->rtrim_nts;
//	argtable[16] = argtable_options->ltrim_nts;
//	argtable[17] = argtable_options->rfilter_nts;
//	argtable[18] = argtable_options->lfilter_nts;
//	argtable[19] = argtable_options->kmers_flag;
//
//	argtable[20] = argtable_options->cg_flag;
//	argtable[21] = argtable_options->k_cg;
//
//	argtable[22] = argtable_options->cpu_num_threads;
//	argtable[23] = argtable_options->batch_size;
//	argtable[24] = argtable_options->batch_list_size;
//
//	argtable[25] = argtable_options->log_level;
//	argtable[26] = argtable_options->log_file;
//	argtable[27] = argtable_options->verbose;
//	argtable[28] = argtable_options->help;
//	argtable[29] = argtable_options->time;
//
//	argtable[30] = end;
//
//	return argtable;
//}
