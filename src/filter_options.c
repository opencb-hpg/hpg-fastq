#include "filter_options.h"

//------------------------------------------------------------------------

void usage_filter_options(filter_options_t *opts);

void **new_argtable_filter_options();
filter_options_t *read_cli_filter_options(void **argtable, filter_options_t *opts);

extern void free_argtable(int num_options, void **argtable);
extern void usage_argtable(char *exec_name, char *command_name, void **argtable);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

filter_options_t *new_filter_options(char *exec_name, char *command_name) {
  filter_options_t *opts = (filter_options_t*) calloc (1, sizeof(filter_options_t));
  
  opts->log_level = 0;
  opts->verbose = 0;
  opts->help = 0;
  opts->num_threads = 2;
  opts->batch_size = 10000;
  opts->in_filename = NULL;
  opts->out_dirname = NULL;
  opts->gff_region_filename = NULL;
  opts->region_list = NULL;

  opts->mapped = 0;
  opts->unmapped = 0;
  opts->unique = 0;
  opts->proper_pairs = 0;
  opts->min_num_errors = -1;
  opts->max_num_errors = -1;
  opts->min_quality = -1;
  opts->max_quality = -1;
  opts->min_length = -1;
  opts->max_length = -1;

  opts->exec_name = strdup(exec_name);
  opts->command_name = strdup(command_name);

  return opts;
}

//------------------------------------------------------------------------

filter_options_t *parse_filter_options(char *exec_name, char *command_name,
				     int argc, char **argv) {
  void **argtable = new_argtable_filter_options();
  
  filter_options_t *opts = new_filter_options(exec_name, command_name);
  if (argc < 2) {
    usage_argtable(exec_name, command_name, argtable);
  } else {  
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*) argtable[0])->count) {
      usage_argtable(exec_name, command_name, argtable);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_FILTER_OPTIONS], exec_name);
      usage_argtable(exec_name, command_name, argtable);
    } else {
      opts = read_cli_filter_options(argtable, opts);
      if (opts->help) {
	usage_argtable(exec_name, command_name, argtable);
      }
    }
  }

  free_argtable(NUM_FILTER_OPTIONS + 1, argtable);

  return opts;
}

//------------------------------------------------------------------------

void free_filter_options(filter_options_t *opts) {
  if (opts == NULL) { return; }
  
  if (opts->in_filename) { free(opts->in_filename); }
  if (opts->out_dirname) { free(opts->out_dirname); }
  if (opts->gff_region_filename) { free(opts->gff_region_filename); }
  if (opts->region_list) { free(opts->region_list); }
  
  if (opts->length_range) { free(opts->length_range); }
  if (opts->quality_range) { free(opts->quality_range); }
  if (opts->num_errors_range) { free(opts->num_errors_range); }

  if (opts->exec_name) { free(opts->exec_name); }
  if (opts->command_name) { free(opts->command_name); }
  
  free(opts);
}

//------------------------------------------------------------------------

void validate_filter_options(filter_options_t *opts) {
  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n\n");
    usage_filter_options(opts);
  }

  if (! exists(opts->out_dirname)) {
    opts->out_dirname = strdup(".");
  }

  int min_value, max_value;
  char min[512], *max;

  if (opts->length_range) {
    strcpy(min, opts->length_range);
    max = strstr(min, "-");
    if (max == NULL) {
      printf("\nError: invalid alignment length range %s: the mininum and maximum value must be separted by a comma, e.g.: 95,105\n\n", opts->length_range);
      usage_filter_options(opts);
    }
    *max = 0;
    opts->min_length = atoi(min);
    opts->max_length = atoi(max + 1);
  }

  if (opts->quality_range) {
    strcpy(min, opts->quality_range);
    max = strstr(min, "-");
    if (max == NULL) {
      printf("\nError: invalid alignment quality range %s: the mininum and maximum value must be separted by a comma, e.g.: 200,245\n\n", opts->quality_range);
      usage_filter_options(opts);
    }
    *max = 0;
    opts->min_quality = atoi(min);
    opts->max_quality = atoi(max + 1);
  }

  if (opts->num_errors_range) {
    strcpy(min, opts->num_errors_range);
    max = strstr(min, "-");
    if (max == NULL) {
      printf("\nError: invalid num. errors range %s: the mininum and maximum value must be separted by a comma, e.g.: 0,3\n\n", opts->num_errors_range);
      usage_filter_options(opts);
    }
    *max = 0;
    opts->min_num_errors = atoi(min);
    opts->max_num_errors = atoi(max + 1);
  }
}

//------------------------------------------------------------------------

void display_filter_options(filter_options_t *opts) {
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Main options\n");
  printf("\tBAM input filename  : %s\n", opts->in_filename);
  printf("\tOutput dirname      : %s\n", opts->out_dirname);
  if (opts->region_list) {
    printf("\tRegions             : %s\n", opts->region_list);
  } else if (opts->gff_region_filename) {
    printf("\tGFF region filename : %s\n", opts->gff_region_filename);
  }
  printf("\n");

  printf("Filters\n");
  if (opts->mapped) printf("\tBy mapped reads\n");
  if (opts->unmapped) printf("\tBy unmapped reads\n");
  if (opts->unique) printf("\tBy unique alignments\n");
  if (opts->proper_pairs) printf("\tBy proper pairs\n");
  if (opts->quality_range) printf("\tBy quality range        : %i - %i\n", opts->min_quality, opts->max_quality);
  if (opts->length_range) printf("\tBy length range          : %i - %i\n", opts->min_length, opts->max_length);
  if (opts->num_errors_range) printf("\tBy num. errors range : %i - %i\n", opts->min_num_errors, opts->max_num_errors);
  printf("\n");

  printf("Report options\n");
  printf("\tLog level: %d\n",  (int) opts->log_level);
  printf("\tVerbose  : %d\n",  (int) opts->verbose);
  printf("\n");

  printf("Architecture options\n");
  printf("\tNum. threads: %d\n",  (int) opts->num_threads);
  printf("\tBatch size  : %d alignments\n",  (int) opts->batch_size);
  printf("=================================================\n");
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

filter_options_t *read_cli_filter_options(void **argtable, filter_options_t *opts) {	
  if (((struct arg_int*)argtable[0])->count) { opts->help = ((struct arg_int*)argtable[0])->count; }
  if (((struct arg_file*)argtable[1])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { opts->out_dirname = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  if (((struct arg_int*)argtable[3])->count) { opts->log_level = *(((struct arg_int*)argtable[3])->ival); }
  if (((struct arg_int*)argtable[4])->count) { opts->verbose = *(((struct arg_int*)argtable[4])->ival); }
  if (((struct arg_int*)argtable[5])->count) { opts->num_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { opts->batch_size = *(((struct arg_int*)argtable[6])->ival); }

  if (((struct arg_int*)argtable[7])->count) { opts->mapped = ((struct arg_int*)argtable[7])->count; }
  if (((struct arg_int*)argtable[8])->count) { opts->unmapped = ((struct arg_int*)argtable[8])->count; }
  if (((struct arg_int*)argtable[9])->count) { opts->unique = ((struct arg_int*)argtable[9])->count; }
  if (((struct arg_int*)argtable[10])->count) { opts->proper_pairs = ((struct arg_int*)argtable[10])->count; }
  
  if (((struct arg_str*)argtable[11])->count) { opts->length_range = strdup(*(((struct arg_str*)argtable[11])->sval)); }
  if (((struct arg_str*)argtable[12])->count) { opts->quality_range = strdup(*(((struct arg_str*)argtable[12])->sval)); }
  if (((struct arg_str*)argtable[13])->count) { opts->num_errors_range = strdup(*(((struct arg_str*)argtable[13])->sval)); }
  
  //  if (((struct arg_file*)argtable[14])->count) { opts->gff_region_filename = strdup(*(((struct arg_file*)argtable[14])->filename)); }
  //  if (((struct arg_str*)argtable[15])->count) { opts->region_list = strdup(*(((struct arg_str*)argtable[15])->sval)); }

  return opts;
}

//--------------------------------------------------------------------

void usage_filter_options(filter_options_t *opts) {
  void **argtable = new_argtable_filter_options();
  usage_argtable(opts->exec_name, opts->command_name, argtable);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void** new_argtable_filter_options() {
  void **argtable = (void**)malloc((NUM_FILTER_OPTIONS + 1) * sizeof(void*));
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  argtable[0] = arg_lit0("h", "help", "Help option");
  argtable[1] = arg_file0("b", "bam-file", NULL, "Input file name (BAM format)");
  argtable[2] = arg_file0("o", "outdir", NULL, "Output file name (BAM format)");
  argtable[3] = arg_int0(NULL, "log-level", NULL, "Log debug level");
  argtable[4] = arg_int0("v", "verbose", NULL, "Verbose");
  argtable[5] = arg_int0(NULL, "num-threads", NULL, "Number of threads");
  argtable[6] = arg_int0(NULL, "batch-size", NULL, "Batch size (in number of alignments)");
  
  argtable[7] = arg_lit0(NULL, "mapped", "Filter by mapped reads");
  argtable[8] = arg_lit0(NULL, "unmapped", "Filter by unmapped reads");
  argtable[9] = arg_lit0(NULL, "unique", "Filter by unique alignments");
  argtable[10] = arg_lit0(NULL, "proper-pairs", "Filter by proper pairs");

  argtable[11] = arg_str0(NULL, "alignment-size", NULL, "Filter by alignment size range (e.g., 95-105)");
  argtable[12] = arg_str0(NULL, "quality", NULL, "Filter by quality range (e.g., 210-250)");
  argtable[13] = arg_str0(NULL, "num-errors", NULL, "Filter by number of errors range (e.g., 0-3)");

  //  argtable[14] = arg_file0(NULL, "gff-refion-file", NULL, "Region file name (GFF format)");
  //  argtable[15] = arg_str0(NULL, "region-list", NULL, "Regions (e.g., 1:3000-3200,4:100-200,...)");

  argtable[NUM_FILTER_OPTIONS] = arg_end(20);
  
  return argtable;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------



/*
int read_config_file(const char *filename, options_t *opts) {
	if (filename == NULL || opts == NULL) {
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

	config_destroy(config);
	free(config);

	return ret_code;
}
*/
