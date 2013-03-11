#include "stats_options.h"

//------------------------------------------------------------------------

void usage_stats_options(stats_options_t *opts);

void **new_argtable_stats_options();
stats_options_t *read_cli_stats_options(void **argtable, stats_options_t *opts);

void free_argtable(int num_options, void **argtable);
void usage_argtable(char *exec_name, char *command_name, void **argtable);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

stats_options_t *new_stats_options(char *exec_name, char *command_name) {
  stats_options_t *opts = (stats_options_t*) calloc (1, sizeof(stats_options_t));
  
  opts->log_level = 0;
  opts->verbose = 0;
  opts->help = 0;
  opts->num_threads = 2;
  opts->batch_size = 10000;
  opts->in_filename = NULL;
  opts->out_dirname = NULL;
  opts->gff_region_filename = NULL;
  opts->region_list = NULL;

  opts->exec_name = strdup(exec_name);
  opts->command_name = strdup(command_name);

  return opts;
}

//------------------------------------------------------------------------

stats_options_t *parse_stats_options(char *exec_name, char *command_name,
				     int argc, char **argv) {
  void **argtable = new_argtable_stats_options();
  
  stats_options_t *opts = new_stats_options(exec_name, command_name);
  if (argc < 2) {
    usage_argtable(exec_name, command_name, argtable);
  } else {  
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*) argtable[0])->count) {
      usage_argtable(exec_name, command_name, argtable);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_STATS_OPTIONS], exec_name);
      usage_argtable(exec_name, command_name, argtable);
    } else {
      opts = read_cli_stats_options(argtable, opts);
      if (opts->help) {
	usage_argtable(exec_name, command_name, argtable);
      }
    }
  }

  free_argtable(NUM_STATS_OPTIONS + 1, argtable);

  return opts;
}

//------------------------------------------------------------------------

void free_stats_options(stats_options_t *opts) {
  if (opts == NULL) { return; }
  
  if (opts->in_filename) { free(opts->in_filename); }
  if (opts->out_dirname) { free(opts->out_dirname); }
  if (opts->gff_region_filename) { free(opts->gff_region_filename); }
  if (opts->region_list) { free(opts->region_list); }
  
  if (opts->exec_name) { free(opts->exec_name); }
  if (opts->command_name) { free(opts->command_name); }
  
  free(opts);
}

//------------------------------------------------------------------------

void validate_stats_options(stats_options_t *opts) {
  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n\n");
    usage_stats_options(opts);
  }

  if (! exists(opts->out_dirname)) {
    opts->out_dirname = strdup(".");
  }
}

//------------------------------------------------------------------------

void display_stats_options(stats_options_t *opts) {
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Main options\n");
  printf("\tFastQ input filename  : %s\n", opts->in_filename);
  printf("\tOutput dirname      : %s\n", opts->out_dirname);
  if (opts->region_list) {
    printf("\tRegions             : %s\n", opts->region_list);
  } else if (opts->gff_region_filename) {
    printf("\tGFF region filename : %s\n", opts->gff_region_filename);
  }
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

stats_options_t *read_cli_stats_options(void **argtable, stats_options_t *opts) {	
  if (((struct arg_int*)argtable[0])->count) { opts->help = ((struct arg_int*)argtable[0])->count; }
  if (((struct arg_file*)argtable[1])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { opts->out_dirname = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  if (((struct arg_int*)argtable[3])->count) { opts->log_level = *(((struct arg_int*)argtable[3])->ival); }
  if (((struct arg_int*)argtable[4])->count) { opts->verbose = *(((struct arg_int*)argtable[4])->ival); }
  if (((struct arg_int*)argtable[5])->count) { opts->num_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { opts->batch_size = *(((struct arg_int*)argtable[6])->ival); }
  /*
  if (((struct arg_file*)argtable[7])->count) { opts->gff_region_filename = strdup(*(((struct arg_file*)argtable[7])->filename)); }
  if (((struct arg_str*)argtable[8])->count) { opts->region_list = strdup(*(((struct arg_str*)argtable[8])->sval)); }
  */
  return opts;
}

//--------------------------------------------------------------------

void usage_stats_options(stats_options_t *opts) {
  void **argtable = new_argtable_stats_options();
  usage_argtable(opts->exec_name, opts->command_name, argtable);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void** new_argtable_stats_options() {
  void **argtable = (void**)malloc((NUM_STATS_OPTIONS + 1) * sizeof(void*));
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  argtable[0] = arg_lit0("h", "help", "Help option");
  argtable[1] = arg_file0("f", "fastq-file", NULL, "Input file name (FastQ format)");
  argtable[2] = arg_file0("o", "outdir", NULL, "Output directory name");
  argtable[3] = arg_int0(NULL, "log-level", NULL, "Log debug level");
  argtable[4] = arg_int0("v", "verbose", NULL, "Verbose");
  argtable[5] = arg_int0(NULL, "num-threads", NULL, "Number of threads");
  argtable[6] = arg_int0(NULL, "batch-size", NULL, "Batch size (in number of alignments)");
  //  argtable[7] = arg_file0(NULL, "gff-refion-file", NULL, "Region file name (GFF format)");
  //  argtable[8] = arg_str0(NULL, "region-list", NULL, "Regions (e.g., 1:3000-3200,4:100-200,...)");
  
  argtable[NUM_STATS_OPTIONS] = arg_end(20);
  
  return argtable;
}

//--------------------------------------------------------------------

void free_argtable(int num_options, void **argtable) {
  if (argtable != NULL) {
    arg_freetable(argtable, num_options);
    free(argtable);
  }
}

//--------------------------------------------------------------------

void usage_argtable(char *exec_name, char *command_name, void **argtable) {
  printf("Usage: %s %s\n", exec_name, command_name);
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
  exit(-1);
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
