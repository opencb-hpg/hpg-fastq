#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "commons.h"
#include "error.h"
#include "log.h"
#include "fastq_file.h"
#include "fastq_read.h"
#include "fastq_batch.h"
#include "qc_batch.h"
#include "qc.h"
#include "file_utils.h"
#include "string_utils.h"
#include "system_utils.h"

#include "cuda_utils.h"

#include "fastq_batch_reader.h"


// global variables for log functions
//
int log_level = DEFAULT_LOG_LEVEL;
bool verbose = true;

// global variables for timing and capacity meausures
//
bool time_on = false;

double read_time = 0.0;
struct timeval t1_read, t2_read;

double gpu_time = 0.0;
struct timeval t1_gpu, t2_gpu;

double result_time = 0.0;
struct timeval t1_result, t2_result;

double reporting_time = 0.0;
struct timeval t1_reporting, t2_reporting;

double total_time = 0.0;
struct timeval t1_total, t2_total;

int number_of_batchs = 0;
double mean_reads_per_batch = 0;
double mean_batch_size = 0;

// FILE* before_reader_fd;
// FILE* reader_fd;

int main(int argc, char **argv) {
  
  // variables to store the value options
  //
  size_t batch_size = 500000; // in bytes
  int num_blocks = 16;
  int num_threads = 256;
  int max_fastq_batch_list_length = 10;
  int kmers_on = 0;
  int base_quality = PHRED33;
  char* fasta_input = (char*)"";
  char* fasta1_input = (char*)"";
  char* fasta2_input = (char*)"";
  char* fastq_input = (char*)"";
  char* fastq1_input = (char*)"";
  char* fastq2_input = (char*)"";
  char* report_directory = (char*)"";
  char* str_quality = (char*)"";
  
  int c;
  int option_index = 0;
	
  // struct defining options and its associated option internal value
  static struct option long_options[] = {
    {"batch-size",     required_argument, 0, 'a'},	 
    {"grid-block-size",           required_argument, 0, 'b'},
    {"fasta",  	      required_argument, 0, 'c'},
    {"fasta1",  	      required_argument, 0, 'd'},
    {"fasta2",  	      required_argument, 0, 'e'},
    {"fastq",  	      required_argument, 0, 'f'},
    {"fq",  	      required_argument, 0, 'f'},
    {"fastq1",  	      required_argument, 0, 'g'},
    {"fq1",  	      required_argument, 0, 'g'},
    {"fastq2",  	      required_argument, 0, 'h'},
    {"fq2",  	      required_argument, 0, 'h'},
    {"report-outdir",  	      required_argument, 0, 'i'},
    {"conf",  	      required_argument, 0, 'j'},
    {"log-level",  	      required_argument, 0, 'k'},
    {"verbose",  	      required_argument, 0, 'l'},
    {"list-length", 	required_argument, 0, 'm'},
    {"phred-quality",  	      required_argument, 0, 'n'},
    {"kmers",  	      no_argument, 0, 'o'},
    {"t",  	      no_argument, 0, '0'},
    {"time",  	      no_argument, 0, '0'},
    {0, 0, 0, 0}
  }; 
  
		
  set_log_level(1);

  int argc_with_file_options;
  char **argv_with_file_options; // = (char **)malloc((argc+2*num_conf_lines) * sizeof(char *));
  char **argv_from_file_options;
	
  for(int i=0; i<argc; i++) {
    if(strcmp(argv[i], "--conf") == 0) {
      
      char str[256];
      strcpy(str, "Reading config file: ");
      strcat(str, argv[i+1]);
      CUDA_LOG_DEBUG(str);
      
      argv_from_file_options = parse_conf_file(argv[i+1]);
      int num_conf_lines = count_lines(argv[i+1]);
      argv_with_file_options = (char **)malloc((argc + 2*num_conf_lines) * sizeof(char *));
      int i = array_concat(argv_with_file_options, argc, (const char**)argv, 2*num_conf_lines, (const char**)argv_from_file_options);
      
      char str1[1024];
      strcpy(str1,"Command line: ");
      argc_with_file_options = argc+2*num_conf_lines;
      
      for(int i=0; i<argc_with_file_options; i++) {
	strcat(str1, argv_with_file_options[i]);
      }
      CUDA_LOG_INFO(str1);
    }
  }
	
// 	if(argv_with_file_options == NULL) {
  argc_with_file_options = argc;
  argv_with_file_options = argv;
// 	}
	
	
  // validation of no argument launch
  //	
  if (argc < 2) {
    printf(QC_USAGE_HELP); 
    exit(0);
  }	
  
  while((c = getopt_long(argc_with_file_options, argv_with_file_options, "", long_options, &option_index)) != -1) {
	
    switch (c) {
    case 'a':
      if (batch_size == 500000) {
	sscanf(optarg, "%li", &batch_size);
      }
      //printf("option --batch-size with value '%li'\n", batch_size);
      break;
      
    case 'b': {
      //printf("option --grid-block-size with value '%s'\n", optarg);
      if (num_blocks == 32 && num_threads == 256) {
	const char delimiters[] = "xX,";
	char *optarg_copy;
	
	optarg_copy = strdupa(optarg);
	sscanf(strtok(optarg_copy, delimiters), "%i", &num_blocks);   // first word => num_blocks
	sscanf(strtok(NULL, delimiters), "%i", &num_threads);   // first word => num_blocks
	//printf("num_blocks: %i, num_threads: %i\n", num_blocks, num_threads);
      }
      break;
    }
    
    case 'c':
      //printf("option --fasta with value '%s'\n", optarg);
      if (strcmp(fasta_input, "") == 0) {
	fasta_input = optarg;
      }
      break;
      
    case 'd':
      //printf("option --fasta1 with value '%s'\n", optarg);
      if (strcmp(fasta1_input, "") == 0) {
	fasta1_input = optarg;
      }
      break;
      
    case 'e':
      //printf("option --fasta2 with value '%s'\n", optarg);
      if (strcmp(fasta2_input, "") == 0) {
	fasta2_input = optarg;
      }
      break;
      
    case 'f':
      //printf("option --fastq or --fq with value '%s'\n", optarg);
      if (strcmp(fastq_input, "") == 0) {
	fastq_input = optarg;
	// 					strcpy(fastq_input, optarg);
      }
      break;
      
    case 'g':
      //printf("option --fastq1 or --fq1 with value '%s'\n", optarg);
      if (strcmp(fastq1_input, "") == 0) {
	fastq1_input = optarg;
      }
      break;
      
    case 'h':
      //printf("option --fastq2 or --fq2 with value '%s'\n", optarg);
      if (strcmp(fastq2_input, "") == 0) {	
	fastq2_input = optarg;
      }
      break;				
      
    case 'i':
      //printf("option --report-dir with value '%s'\n", optarg);
      //if(strcmp(report_directory, "") == 0) {
      report_directory = optarg;
      //}
      break;				
      
    case 'k':
      //printf("option --log-level with value '%s'\n", optarg);
      if (log_level == -1) {
	set_log_level(atoi(optarg));
      }
      break;	
      
      
    case 'l':
      //printf("option --verbose with value '%s'\n", optarg);
      log_set_verbose(strcmp(optarg, "true") == 0);
      break;	
      
    case 'm':
      //printf("option --list-length with value '%s'\n", optarg);
      if (max_fastq_batch_list_length == 10) {
	sscanf(optarg, "%i", &max_fastq_batch_list_length);
      }
      break;	      

    case 'n':
      //printf("option --quality-score with value '%s'\n", optarg);
      if (strcmp(str_quality, "") == 0) {
	str_quality = optarg;
	
	if (strcmp(str_quality, "33") == 0) {
	  base_quality = PHRED33;
	} else if (strcmp(str_quality, "64") == 0) {
	  base_quality = PHRED64;
	} else if (strcmp(str_quality, "sanger") == 0) {
	  base_quality = PHRED33;	
	} else if (strcmp(str_quality, "solexa") == 0) {
	  base_quality = PHRED64;	  
	} else {
	  CUDA_LOG_WARN("Incorrect quality scale (33 or 64). Assuming 33.\n");
	}
      }
      break;	
      
    case 'o':
      //printf("option --kmers selected, kmers calculation enabled\n", optarg);
      kmers_on = 1;
      break;	      
      
    case '0':
      //printf("option --time selected, timing enabled\n");
      time_on = true;
      break;	      
      
    case ':':       /* option without mandatory operand */
      fprintf(stderr, "Option -%c requires an operand\n", optopt);
      // 				errflg++;
      break;
      
    case '?':
      printf(QC_USAGE_HELP); 
      break;
      
    default:
      printf(QC_USAGE_HELP); 
      abort ();
    }
  }
    
  // single-end and paired-end options are exclusive (fasta)
  //
  if ((strcmp(fasta_input, "") != 0 && strcmp(fasta1_input, "") != 0) || (strcmp(fasta_input, "") != 0 && strcmp(fasta2_input, "") != 0)) {
    printf("single-end and paired-end options are exclusive, use --fasta OR --fasta1/--fasta2 options, not both\n");    
    printf(QC_USAGE_HELP); 
    exit(0);    
  }	
  
  // single-end and paired-end options are exclusive (fastq)
  if ((strcmp(fastq_input, "") != 0 && strcmp(fastq1_input, "") != 0) || (strcmp(fastq_input, "") != 0 && strcmp(fastq2_input, "") != 0)) {
    printf("single-end and paired-end options are exclusive, use --fastq OR --fastq1/--fastq2 options, not both\n");    
    printf(QC_USAGE_HELP); 
    exit(0);    
  }
  
  // batch-size > 8192 for security (at least 1 read must be contained in the batch)
  if (batch_size < 8192) {
    batch_size = 8192;
    printf("the value of N in --batch-size N must be at least 8192 \n");        
  }

  // fasta and fastq inputs are exclusive
  char *fasta_options = (char*)malloc(strlen(fasta_input)+strlen(fasta1_input)+strlen(fasta2_input));
  char *fastq_options = (char*)malloc(strlen(fastq_input)+strlen(fastq1_input)+strlen(fastq2_input));
  
  strcat(fasta_options, fasta_input);
  strcat(fasta_options, fasta1_input);
  strcat(fasta_options, fasta2_input);

  strcat(fastq_options, fastq_input);
  strcat(fastq_options, fastq1_input);
  strcat(fastq_options, fastq2_input);
  
  //printf("fasta_options: %s\n", fasta_options);
  //printf("fastq_options: %s\n", fastq_options);
  
  if (strcmp(fasta_options, "") != 0 && strcmp(fastq_options, "") != 0) {
	  printf("fasta and fastq inputs are exclusive, select only a type of file\n");    
	  printf(QC_USAGE_HELP); 
	  exit(0);    
  }	
  
  // validation of the mandatory options. getopt_long does not provide this validation
  //
  // at least one fasta/fastq input is mandatory
  //
  if (strcmp(fasta_input, "") == 0 && strcmp(fastq_input, "") == 0 && strcmp(fastq1_input, "") == 0 && strcmp(fastq2_input, "") == 0) {
	  printf("at least one fasta/fastq file input is is mandatory\n");    
	  printf(QC_USAGE_HELP); 
	  exit(0);    
  }

  // both pair ends must be informed
  if ((strcmp(fastq1_input, "") != 0 && strcmp(fastq2_input, "") == 0) || (strcmp(fastq1_input, "") == 0 && strcmp(fastq2_input, "") != 0)) {
	  printf("both pair ends files are mandatory, use both --fastq1 and --fastq2 options\n");    
	  printf(QC_USAGE_HELP); 
	  exit(0);    
  }

  // directory of qc report output is mandatory
  if (strcmp(report_directory, "") == 0) {
	  printf("--report-dir is mandatory\n");    
	  printf(QC_USAGE_HELP); 
	  exit(0);
  }  
  
  // free memory must be enough to launch quality control
  if (get_estimated_memory_needed(QC, batch_size, max_fastq_batch_list_length) > get_free_memory()) {
	  set_log_level(1);
	  verbose = true;
	  CUDA_LOG_ABORT("There is no enough free memory to launch the program. Please reduce batch size or list length.\n");
  }
	
  // print any remaining command line arguments that are not options
  //
  if (optind < argc) {
	  printf("no valid options: ");
	  while (optind < argc) printf("%s ", argv[optind++]);
	  printf(QC_USAGE_HELP);
  }

  // resources must be free
  free(fasta_options);
  free(fastq_options);
  
  if (time_on) { start_timer(t1_total); }

  //file descriptors to obtain reads during qc intermediate phases
  //
//   before_reader_fd = fopen("/tmp/before_reader.fastq", "w");
//   reader_fd = fopen("/tmp/reader.fastq", "w");


  num_blocks = (batch_size / 2 / num_threads) + 1;
  
  if (strcmp(fastq_input, "") != 0) {
    kernel_qc_fastq_single_end(batch_size, max_fastq_batch_list_length, num_blocks, num_threads, base_quality, kmers_on, fastq_input, report_directory);
  } else if (strcmp(fastq1_input, "") != 0 && strcmp(fastq2_input, "") != 0) {
    kernel_qc_fastq_paired_end(batch_size, max_fastq_batch_list_length, num_blocks, num_threads, base_quality, kmers_on, fastq1_input, fastq2_input, report_directory);
  } else {
    printf("Missing input files\n");
    printf(QC_USAGE_HELP);
  }
  
//   fclose(before_reader_fd);
//   fclose(reader_fd);
  
  total_time = 0;
  if (time_on) { stop_timer(t1_total, t2_total, total_time); }

  if (time_on) {
   
    printf("\n");  
    printf("number of batches     : \t%10i\n\n", number_of_batchs);
    printf("mean reads per batch  : \t%10.2f\n", mean_reads_per_batch / number_of_batchs); 
    printf("mean batch size (KB)  : \t%10.2f\n\n", mean_batch_size / number_of_batchs / 1024); 
    
    printf("total time           (s): \t%10.5f\n", 0.000001 * total_time);  
    printf("\n");  
    printf("total read time      (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * read_time, 0.000001 * read_time / number_of_batchs);  
    printf("total gpu time       (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * gpu_time, 0.000001 * gpu_time / number_of_batchs);  
    printf("total result time    (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * result_time, 0.000001 * result_time / number_of_batchs);  
    printf("total reporting time (s): \t%10.5f\n", 0.000001 * reporting_time);  
  }
    
  return 1;
}