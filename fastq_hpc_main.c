#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "chaos_game.h"
#include "commons.h"
#include "cuda_commons.h"
#include "log.h"
#include "file_utils.h"
#include "string_utils.h"
#include "system_utils.h"

#define DEFAULT_MIN_READ_LENGTH    50
#define DEFAULT_MAX_READ_LENGTH    200
#define DEFAULT_MIN_READ_QUALITY    20
#define DEFAULT_MAX_READ_QUALITY    60

#define DEFAULT_GPU_NUM_BLOCKS     16
#define DEFAULT_GPU_NUM_THREADS    512
#define DEFAULT_GPU_NUM_DEVICES    0
#define DEFAULT_CPU_NUM_THREADS    2
#define DEFAULT_BATCH_SIZE_MB     64
#define DEFAULT_BATCH_LIST_SIZE    4

#define DEFAULT_K_IN_CHAOS_GAME    10

#define FASTQ_HPC_TOOLS_USAGE "USAGE: fastq-hpc-tools [--prepro] [--filter] [--qc] --outdir --min-read-length --max-read-length [--last-nts]  \
 [--first-nts] [--begin-quality-nt] [--end-quality-nt] [--min-quality] [--max-quality] [--max-nts-mismatch] [--max-N-per-read] [--batch-size] \
 [--fastq|--fq]|[--fastq1|--fq1 --fastq2|--fq2] [--batch-list-size] [--phred-quality] [--kmers] [--rtrim] [--ltrim] [--t | --time] \
 [--cg -k --gs-filename]"

/* **********************************************
 *    		Global variables  		*
 * **********************************************/

int time_flag = 0;

double read_time = 0.0;
struct timeval t1_read, t2_read;

double gpu_time = 0.0;
struct timeval t1_gpu, t2_gpu;

double cpu_time = 0.0;
struct timeval t1_cpu, t2_cpu;

double kmers_time = 0.0;
struct timeval t1_kmers, t2_kmers;

double result_time = 0.0;
struct timeval t1_result, t2_result;

double write_time = 0.0;
struct timeval t1_write, t2_write;

double chaos_game_time = 0.0;
struct timeval t1_chaos_game, t2_chaos_game;

double reporting_time = 0.0;
struct timeval t1_reporting, t2_reporting;

double total_time = 0.0;
struct timeval t1_total, t2_total;

int number_of_batchs = 0;
double mean_reads_per_batch = 0;
double mean_batch_size = 0;

int main(int argc, char **argv) {
    // setting global variables for logger
    log_level = LOG_DEFAULT_LEVEL;
    verbose = 1;
    log_filename = NULL;

    int qc_flag = 0;
    int kmers_flag = 0;
    int cg_flag = 0; //Chaos Game flag

    int filter_flag = 0;
    int prepro_flag = 0;
    int min_read_length =  DEFAULT_MIN_READ_LENGTH;  // 50
    int max_read_length =  DEFAULT_MAX_READ_LENGTH;  // 200
    int min_quality =  DEFAULT_MIN_READ_QUALITY;  // 20
    int max_quality =  DEFAULT_MAX_READ_QUALITY;  // 60
    int max_nts_mismatch = 3;
    int max_n_per_read = 0;
    int begin_quality_nt = 0;
    int end_quality_nt = 1024;
    int base_quality = PHRED33;
    int rtrim_nts = 0;
    int ltrim_nts = 0;
    int rfilter_nts = 0;
    int lfilter_nts = 0;

    int k_cg = DEFAULT_K_IN_CHAOS_GAME;    //7

    // variables to store the value options
    int gpu_num_blocks =  DEFAULT_GPU_NUM_BLOCKS;   // 16
    int gpu_num_threads =  DEFAULT_GPU_NUM_THREADS;  // 512
    int gpu_num_devices =  DEFAULT_GPU_NUM_DEVICES;  // -1
    int cpu_num_threads =  DEFAULT_CPU_NUM_THREADS;  // 2
    size_t batch_size =  DEFAULT_BATCH_SIZE_MB * 1000000; // 64MB
    int batch_list_size =  DEFAULT_BATCH_LIST_SIZE;  // 4

    // conf file is parsed below
    char* fastq_input = NULL;
    char* fastq1_input = NULL;
    char* fastq2_input = NULL;
    char* genomic_signature_input = NULL;
    char* output_directory = NULL;

    int c;
    int option_index = 0;

    // struct defining options and its associated option internal value
    static struct option long_options[] = {
        /* QC parameters */
        {"qc",    no_argument, 0, 'a'},
        {"quality-control",   no_argument, 0, 'a'},
        {"kmers",    no_argument, 0, 'b'},

        /* Filter and Preprocessing parameters */
        {"filter",    no_argument, 0, 'c'},
        {"prep",    no_argument, 0, 'd'},
        {"preprocessing",   no_argument, 0, 'd'},
        {"min-read-length",   required_argument, 0, 'e'},
        {"max-read-length",   required_argument, 0, 'f'},
        {"min-quality",   required_argument, 0, 'g'},
        {"max-quality",   required_argument, 0, 'h'},
        {"max-nts-out-quality",  required_argument, 0, 'i'},
        {"max-N-per-read",   required_argument, 0, 'j'},
        {"max-n-per-read",   required_argument, 0, 'j'},
        {"start-quality-nt",   required_argument, 0, 'k'}, // By default nt 0
        {"end-quality-nt",   required_argument, 0, 'l'},
        {"rtrim-nts",    required_argument, 0, 'm'},
        {"ltrim-nts",    required_argument, 0, 'n'},
        {"rfilter-nts",   required_argument, 0, 'o'},
        {"lfilter-nts",   required_argument, 0, 'p'},

        /* Commons parameters */
        {"phred-quality",   required_argument, 0, 'q'},

        /* HPC parameters  */
        {"gpu-num-blocks",   required_argument, 0, 'r'}, // num-blocks: 16 (dafault 16)
        {"gpu-num-threads",   required_argument, 0, 's'}, // num-threads: 512
        {"gpu-num-devices",   required_argument, 0, 't'}, // num-devices: 1 (NOTE: 0 for all of them)
        {"cpu-num-threads",   required_argument, 0, 'u'}, // num-threads: 12 (default: 2, all possible when 0)
        {"batch-size",   required_argument, 0, 'v'}, // batch-size (MB): default 256 (~256MB) (NOTE: optimize for GPUs)
        {"batch-list-size",    required_argument, 0, 'w'}, // batch-list-size: number of batches in the list (default: 4)

        /* IO parameters */
        {"fq",    required_argument, 0, 'x'},
        {"fastq",    required_argument, 0, 'x'},
        {"fq1",     required_argument, 0, 'y'},
        {"fastq1",    required_argument, 0, 'y'},
        {"fq2",     required_argument, 0, 'z'},
        {"fastq2",    required_argument, 0, 'z'},
        {"o",     required_argument, 0, 'A'},
        {"outdir",    required_argument, 0, 'A'},
        {"conf",    required_argument, 0, 'B'},

        /* LOG parameters  */
        {"log-level",    required_argument, 0, 'C'}, // levels form 1 to 5 (DEBUG to FATAL)
        {"log-file",    required_argument, 0, 'D'},
        {"v",     required_argument, 0, 'E'}, // verbose: if False no 'console' output
        {"verbose",    required_argument, 0, 'E'},
        {"t",     no_argument, 0, 'F'},
        {"time",    no_argument, 0, 'F'},

        /* GENOMIC SIGNATURE (chaos game) parameters  */
        {"cg",    no_argument, 0, 'G'},
        {"chaos-game",   no_argument, 0, 'G'},
        {"k",     required_argument, 0, 'H'},
        {"gs-filename",   required_argument, 0, 'I'},

        {0, 0, 0, 0}
    };

    LOG_LEVEL(LOG_INFO_LEVEL);

    // parsing confing file if passed as argument (--conf <file>)
    int argc_with_file_options;
    char **argv_with_file_options = NULL;
    char **argv_from_file_options = NULL;

    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--conf") == 0) {
            char* str = (char*) calloc(256, sizeof(char));
            strcpy(str, "Reading config file: ");
            strcat(str, argv[i+1]);
            LOG_DEBUG(str);

            argv_from_file_options = parse_conf_file(argv[i+1]);
            int num_conf_lines = count_lines(argv[i+1]);
            argv_with_file_options = (char **)malloc((argc + 2 * num_conf_lines) * sizeof(char *));
            int i = array_concat(argv_with_file_options, argc, (const char**)argv, 2 * num_conf_lines, (const char**)argv_from_file_options);

            char str1[1024];
            strcpy(str1, "Command line: ");
            argc_with_file_options = argc + 2 * num_conf_lines;

            for (int j = 0; j < argc_with_file_options; j++) {
                strcat(str1, argv_with_file_options[j]);
            }

            LOG_INFO(str1);

            free(str);
        }
    }

    argc_with_file_options = argc;
    argv_with_file_options = argv;

    // validation of no argument launch
    if (argc < 2) {
        printf(FASTQ_HPC_TOOLS_USAGE);
        exit(0);
    }

    while ((c = getopt_long(argc_with_file_options, argv_with_file_options, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 'a':
                qc_flag = 1;
                break;

            case 'b':
                kmers_flag = 1;
                break;

            case 'c':
                filter_flag = 1;
                break;

            case 'd':
                prepro_flag = 1;
                break;

            case 'e':
                //printf("option --min-read-length with value '%s'\n", optarg);
                if (min_read_length == DEFAULT_MIN_READ_LENGTH) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &min_read_length);
                        if (min_read_length < 20) {
                            min_read_length = 20;
                            LOG_WARN("--min-read-length is not a valid number, assuming 20\n");
                        }
                    } else {
                        LOG_FATAL("--min-read-length is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'f':
                //printf("option --max-read-length with value '%s'\n", optarg);
                if (max_read_length == DEFAULT_MAX_READ_LENGTH) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &max_read_length);
                        if (max_read_length <= 0) {
                            LOG_FATAL("--max-read-length is not a valid number, aborting execution\n");
                        }
                    } else {
                        LOG_FATAL("--max-read-length is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'g':
                //printf("option --min-quality with value '%s'\n", optarg);
                if (min_quality == DEFAULT_MIN_READ_QUALITY) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &min_quality);
                        min_quality = max(min_quality, MIN_QUALITY_VALUE);
                    } else {
                        LOG_FATAL("--min-quality is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'h':
                //printf("option --max-quality with value '%s'\n", optarg);
                if (max_quality == DEFAULT_MAX_READ_QUALITY) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &max_quality);
                        max_quality = min(max_quality, MAX_QUALITY_VALUE);
                    } else {
                        LOG_FATAL("--max-quality is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'i':
                //printf("option --max-nts-mismatch with value '%s'\n", optarg);
                if (max_nts_mismatch == 3) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &max_nts_mismatch);
                    } else {
                        LOG_FATAL("--max-nts-mismatch is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'j':
                //printf("option --max-n-per-read with value '%s'\n", optarg);
                if (max_n_per_read == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &max_n_per_read);
                    } else {
                        LOG_FATAL("--max-n-per-read is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'k':
                //printf("option --begin-quality-nt with value '%s'\n", optarg);
                if (begin_quality_nt == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &begin_quality_nt);
                    } else {
                        LOG_FATAL("--begin-quality-nt is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'l':
                //printf("option --end-quality-nt with value '%s'\n", optarg);
                if (end_quality_nt == 1024) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &end_quality_nt);
                    } else {
                        LOG_FATAL("--end-quality-nt is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'm':
                //printf("option --rtrim-nts with value '%s'\n", optarg);
                if (rtrim_nts == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &rtrim_nts);
                    } else {
                        LOG_FATAL("--rtrim-nts is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'n':
                //printf("option --ltrim-nts with value '%s'\n", optarg);
                if (ltrim_nts == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &ltrim_nts);
                    } else {
                        LOG_FATAL("--ltrim-nts is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'o':
                //printf("option --rfilter-nts with value '%s'\n", optarg);
                if (rfilter_nts == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &rfilter_nts);
                    } else {
                        LOG_FATAL("--rfilter-nts is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'p':
                //printf("option --lfilter-nts with value '%s'\n", optarg);
                if (lfilter_nts == 0) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &lfilter_nts);
                    } else {
                        LOG_FATAL("--lfilter-nts is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'q':
                //printf("option --phred-quality with value '%s'\n", optarg);
                if (base_quality == PHRED33) {

                    if (strcmp(optarg, "33") == 0) {
                        base_quality = PHRED33;
                    } else if (strcmp(optarg, "64") == 0) {
                        base_quality = PHRED64;
                    } else if (strcmp(optarg, "sanger") == 0) {
                        base_quality = PHRED33;
                    } else if (strcmp(optarg, "solexa") == 0) {
                        base_quality = PHRED64;
                    } else {
                        LOG_WARN("Incorrect quality scale (33 or 64). Assuming 33.\n");
                    }
                }
                break;

            /* PARSING HPC PARAMETERS */
            case 'r':
                if (gpu_num_blocks == DEFAULT_GPU_NUM_BLOCKS) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &gpu_num_blocks);
                        if (gpu_num_blocks < 8) {
                            gpu_num_blocks = 8;
                            LOG_WARN("--gpu-num-blocks is not a valid number, assuming 8\n");
                        }
                    } else {
                        LOG_FATAL("--gpu-num-blocks is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 's':
                if (gpu_num_threads == DEFAULT_GPU_NUM_THREADS) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &gpu_num_threads);
                        if (gpu_num_threads < 32) {
                            gpu_num_threads = 32;
                            LOG_WARN("--gpu-num-threads is not a valid number, assuming 32\n");
                        }
                    } else {
                        LOG_FATAL("--gpu-num-threads is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 't':
                //printf("option --grid-block-size with value '%s'\n", optarg);
                if (gpu_num_devices == DEFAULT_GPU_NUM_DEVICES) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &gpu_num_devices);
                    } else {
                        LOG_FATAL("--gpu-num-devices is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'u':
                // if gets a different value from conf file we do not enter
                if (cpu_num_threads == DEFAULT_CPU_NUM_THREADS) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &cpu_num_threads);
                    } else {
                        LOG_FATAL("--cpu-num-threads is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'v':
                // if gets a different value from conf file we do not enter
                if (batch_size == DEFAULT_BATCH_SIZE_MB * 1000000) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &batch_size);

                        // batch-size > 64MB
                        if (batch_size < 64) {
                            batch_size = 64000000;
                            LOG_WARN("the value of N in --batch-size N must be at least 64MB \n");
                        }
                    } else {
                        batch_size = 64;
                        LOG_FATAL("--batch-size is not a valid number, aborting execution\n");
                    }
                }
                break;

            case 'w':
                if (batch_list_size == DEFAULT_BATCH_LIST_SIZE) {
                    if (is_numeric(optarg) != 0) {
                        sscanf(optarg, "%i", &batch_list_size);

                        if (batch_list_size < 2) {
                            batch_list_size = 2;
                            LOG_WARN("the value of N in --batch-list-size N must be at least 2 \n");
                        }
                    } else {
                        LOG_FATAL("--batch-list-size is not a valid number, aborting execution, assuming default value 4\n");
                    }
                }
                break;

            /* PARSING IO PARAMETERS */
            case 'x':
                //printf("option --fastq or --fq with value '%s'\n", optarg);
                if (fastq_input == NULL) {
                    fastq_input = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(fastq_input, optarg);
                }
                break;

            case 'y':
                //printf("option --fastq1 or --fq1 with value '%s'\n", optarg);
                if (fastq1_input == NULL) {
                    fastq1_input = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(fastq1_input, optarg);
                }
                break;

            case 'z':
                //printf("option --fastq2 or --fq2 with value '%s'\n", optarg);
                if (fastq2_input == NULL) {
                    fastq2_input = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(fastq2_input, optarg);
                }
                break;

            case 'A':
                //printf("option --outdir with value '%s'\n", optarg);
                if (output_directory == NULL) {
                    output_directory = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(output_directory, optarg);
                }
                break;

            /* PARSING LOG PARAMETERS */
            case 'C':
                //printf("option --log-level with value '%s'\n", optarg);
                if (is_numeric(optarg) != 0) {
                    sscanf(optarg, "%i", &log_level);
                    LOG_LEVEL(log_level);
                } else {
                    LOG_WARN("--log-level is not a valid number, assuming default level ERROR\n");
                }
                break;

            case 'D':
                //printf("option --log-file with value '%s'\n", optarg);
                if (log_filename == NULL) {
                    log_filename = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(log_filename, optarg);
                }
                break;

            case 'E':
                //printf("option --verbose or --v with value '%s'\n", optarg);
                if (strcmp(optarg, "true") == 0) {
                    LOG_VERBOSE(1);
                } else if (strcmp(optarg, "false") == 0) {
                    LOG_VERBOSE(0);
                } else {
                    LOG_WARN("--verbose parameter must be true or false, assuming default value false");
                }
                break;

            case 'F':
                //printf("option --time selected. timing enabled\n");
                time_flag = 1;
                break;

            case 'G':
                //printf("option --cg or --chaos-game selected. genomic signature (chaos game) calculation enabled\n");
                cg_flag = 1;
                break;

            case 'H':
                //printf("option --k for chaos game with value %i\n", optarg);
                if (is_numeric(optarg) != 0) {
                    sscanf(optarg, "%i", &k_cg);
                } else {
                    LOG_WARN("--k is not a valid number, assuming default value 7\n");
                }
                break;

            case 'I':
                //printf("option --genomic_signature_input with value '%s'\n", optarg);
                if (genomic_signature_input == NULL) {
                    genomic_signature_input = (char*) calloc(strlen(optarg) + 1, sizeof(char));
                    strcpy(genomic_signature_input, optarg);
                }
                break;

            case ':':       /* option without mandatory operand */
                fprintf(stderr, "Option -%c requires an operand\n", optopt);
                break;

            case '?':
                printf(FASTQ_HPC_TOOLS_USAGE);
                break;

            default:
                printf(FASTQ_HPC_TOOLS_USAGE);
                abort();
        }
    }

    // de-normalization of quality scale
    min_quality += base_quality;
    max_quality += base_quality;

    // fasta and fastq inputs are exclusive
    int fastq_options_length = 3;
    if (fastq_input != NULL) fastq_options_length += strlen(fastq_input);
    if (fastq1_input != NULL) fastq_options_length += strlen(fastq1_input);
    if (fastq2_input != NULL) fastq_options_length += strlen(fastq2_input);

    char *fastq_options = (char*) calloc(fastq_options_length, sizeof(char));

    if (fastq_input != NULL) strcat(fastq_options, fastq_input);
    if (fastq1_input != NULL) strcat(fastq_options, fastq1_input);
    if (fastq2_input != NULL) strcat(fastq_options, fastq2_input);

    // single-end and paired-end options are exclusive (fastq)
    if (fastq_input != NULL && fastq1_input != NULL || fastq_input != NULL && fastq2_input != NULL) {
        printf(FASTQ_HPC_TOOLS_USAGE);
        LOG_FATAL("single-end and paired-end options are exclusive, use --fastq OR --fastq1/--fastq2 options, not both\n");
    }

    // validation of the mandatory options. getopt_long does not provide this validation
    // at least one fasta/fastq input is mandatory
    if (fastq_input != NULL) {
        char log_message[100];
        sprintf(log_message, "fastq file to process: %s\n", fastq_input);
        LOG_DEBUG(log_message);
    } else {
        // both pair ends must be informed
        if (fastq1_input == NULL  || fastq2_input == NULL) {
            printf("Both pair ends files are mandatory, use both --fastq1 and --fastq2 options\n");
            printf(FASTQ_HPC_TOOLS_USAGE);
            exit(0);
        }
    }

    // directory of files output is mandatory
    if (output_directory == NULL) {
        LOG_ERROR("--outdir is mandatory\n");
        printf(FASTQ_HPC_TOOLS_USAGE);
        exit(0);
    }


    // if not flag is passed then we assume QC
    if ((prepro_flag == 0) && (filter_flag == 0) && (qc_flag == 0)) {
        qc_flag = 1;
    }

    // aplying heuristic values if default values have not been modified
    if (batch_size == 1000000 * DEFAULT_BATCH_SIZE_MB) {
        if ((prepro_flag != 0) || (filter_flag != 0)) {
            batch_size = get_optimal_batch_size(FASTQ_PREPRO, batch_list_size);
        } else {
            batch_size = get_optimal_batch_size(FASTQ_QC, batch_list_size);
        }
    }

    if (cpu_num_threads == DEFAULT_CPU_NUM_THREADS) {
        cpu_num_threads = get_optimal_cpu_num_threads();
    }

    if (gpu_num_threads == DEFAULT_GPU_NUM_THREADS) {
        gpu_num_threads = get_optimal_gpu_num_threads();
    }

    // number of rtrim nucleotides cannot be more than half the minumum read length
    if (rtrim_nts > (min_read_length / 4)) {
        LOG_ERROR("--rtrim-nts must be at most 1/4 the value of min_read_length\n");
        exit(0);
    }

    // number of ltrim nucleotides cannot be more than half the minumum read length
    if (ltrim_nts > (min_read_length / 4)) {
        LOG_ERROR("--ltrim-nts must be at most 1/4 the value of min_read_length\n");
        exit(0);
    }

    // number of rfilter nucleotides cannot be more than half the minumum read length
    if (rfilter_nts > (min_read_length / 4)) {
        LOG_ERROR("--rfilter-nts must be at most 1/4 the value of min_read_length\n");
        exit(0);
    }

    // number of lfilter nucleotides cannot be more than half the minumum read length
    if (lfilter_nts > (min_read_length / 4)) {
        LOG_ERROR("--lfilter-nts must be at most 1/4 the value of min_read_length\n");
        exit(0);
    }

    // print any remaining command line arguments that are not options
    if (optind < argc) {
        LOG_WARN("no valid options: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        printf("\n");
        printf(FASTQ_HPC_TOOLS_USAGE);
    }

    // start timer for total time
    if (time_flag) {
        start_timer(t1_total);
    }

    gpu_num_blocks = (batch_size / 2 / gpu_num_threads) + 1;

    if (fastq_input != NULL) {
        kernel_prepro_fastq_single_end(batch_size, batch_list_size, gpu_num_blocks, gpu_num_threads, cpu_num_threads, fastq_input, output_directory, min_quality, max_quality, base_quality, begin_quality_nt, end_quality_nt, max_nts_mismatch, max_n_per_read, min_read_length, max_read_length, rtrim_nts, ltrim_nts, rfilter_nts, lfilter_nts, prepro_flag, filter_flag, qc_flag, kmers_flag, cg_flag, k_cg, genomic_signature_input);
    } else if (fastq1_input != NULL && fastq2_input != NULL) {
        kernel_prepro_fastq_paired_end(batch_size, batch_list_size, gpu_num_blocks, gpu_num_threads, cpu_num_threads, fastq1_input, fastq2_input, output_directory, min_quality, max_quality, base_quality, begin_quality_nt, end_quality_nt, max_nts_mismatch, max_n_per_read, min_read_length, max_read_length, rtrim_nts, ltrim_nts, rfilter_nts, lfilter_nts, prepro_flag, filter_flag, qc_flag, kmers_flag, cg_flag, k_cg, genomic_signature_input);
    } else {
        printf("Missing input files\n");
        printf(FASTQ_HPC_TOOLS_USAGE);
    }

    // stop timer for total time
    total_time = 0;
    if (time_flag) {
        stop_timer(t1_total, t2_total, total_time);
    }

    // print timer measures
    if (time_flag) {
        printf("\n");
        printf("number of batches     : \t%10i\n\n", number_of_batchs);
        printf("mean reads per batch  : \t%10.2f\n", mean_reads_per_batch / number_of_batchs);
        printf("mean batch size (KB)  : \t%10.2f\n\n", (mean_batch_size / number_of_batchs) / 1024);

        printf("total time            (s): \t%10.5f\n", 0.000001 * total_time);
        printf("\n");
        printf("total read time       (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * read_time, 0.000001 * read_time / number_of_batchs);
        printf("total gpu time        (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * gpu_time, 0.000001 * gpu_time / number_of_batchs);
        printf("total cpu time        (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * cpu_time, 0.000001 * cpu_time / number_of_batchs);
        printf("total kmers time      (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * kmers_time, 0.000001 * kmers_time / number_of_batchs);
        printf("total chaos game time (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * chaos_game_time, 0.000001 * chaos_game_time / number_of_batchs);
        printf("total result time     (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * result_time, 0.000001 * result_time / number_of_batchs);

        if ((prepro_flag) || (filter_flag)) {
            printf("total write time      (s): \t%10.5f\t\tper batch: %10.5f\n", 0.000001 * write_time, 0.000001 * write_time / number_of_batchs);
        }

        if (qc_flag) {
            printf("total reporting time  (s): \t%10.5f\n", 0.000001 * reporting_time);
        }
    }

    //free memory for variables
    if (log_filename != NULL) free(log_filename);
    if (argv_from_file_options != NULL) free(argv_from_file_options);
    if (fastq_input != NULL) free(fastq_input);
    if (fastq1_input != NULL) free(fastq1_input);
    if (fastq2_input != NULL) free(fastq2_input);
    if (output_directory != NULL) free(output_directory);
    if (genomic_signature_input != NULL) free(genomic_signature_input);
    if (fastq_options != NULL) free(fastq_options);

    return 1;
}
