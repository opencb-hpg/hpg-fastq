
#ifndef PREPRO_CU
#define PREPRO_CU

#include "chaos_game.h"
#include "fastq_batch_reader.h"
#include "fastq_read.h"
#include "qc_batch.h"
#include "prepro.h"
#include "prepro_batch.h"
#include "prepro_kernel_omp.h"


#define BLOCK_SIZE  16

/* **********************************************
 *    		Private functions  		*
 * *********************************************/

void* qc_calc_server(void* param_p);
void* results_server(void* param_p);
void* writer_single_end_server(void* param_p);
void* writer_paired_end_server(void* param_p);

/* **********************************************
 *    		Global variables   		*
 * *********************************************/

// the following global variables for reads and results batchs are critical structures if used within threads
// they must be accessed EXCLUSIVELY when modifying
list_t qc_batch_list;
status_batch_list_t status_batch_list;

int gpus_thread_alive = 1;
int results_thread_alive = 1;

pthread_mutex_t gpus_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t results_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;

/* **************************************************************
 *    		Private functions implementations  		*
 * **************************************************************/

/* **********************************************
 *    		QC calculator server   		*
 * *********************************************/

void* qc_calc_server(void* params_p) {
    LOG_DEBUG("Thread-GPU: START\n");
    qc_calc_server_input_t* input_p = (qc_calc_server_input_t*) params_p;

    int cpu_num_threads = input_p->cpu_num_threads;
    fastq_batch_list_item_t* fastq_batch_list_item_p = NULL;
    fastq_batch_list_t* batch_list_p = input_p->batch_list_p;
    qc_batch_t* qc_batch_p = NULL;
    list_item_t* item_p = NULL;
    int num_total_threads = input_p->num_blocks * input_p->num_threads;

    // variables for store output results in both CPU and GPU
    qc_read_t* gpu_result_p = NULL;
    qc_kmers_t* gpu_kmers_valid_p = NULL;
    qc_kmers_t* gpu_kmers_invalid_p = NULL;

    char* d_seq_p;
    char* d_quality_p;
    int* d_data_indices_p;
    qc_read_t* d_gpu_result_p;

    int reads_alive;

    reads_alive = fastq_batch_list_get_producers(batch_list_p);

    fastq_batch_list_item_p = fastq_batch_list_remove(batch_list_p);

    while (reads_alive > 0 || fastq_batch_list_item_p != NULL) {
        if (fastq_batch_list_item_p == NULL) {
            // Delay for a bit
            sched_yield();
            usleep(1000);
        } else {
            char log_message[100];
            sprintf(log_message, "Thread-GPU: processing for batch %i....\n", fastq_batch_list_item_p->id);
            LOG_DEBUG(log_message);

            // allocation memory for output results
            gpu_result_p = (qc_read_t*) calloc(fastq_batch_list_item_p->batch_p->num_reads, sizeof(qc_read_t));
            if (input_p->kmers_on) {
                gpu_kmers_valid_p = (qc_kmers_t*) malloc(KMERS_COMBINATIONS * sizeof(qc_kmers_t));
                gpu_kmers_invalid_p = (qc_kmers_t*) malloc(KMERS_COMBINATIONS * sizeof(qc_kmers_t));
                memset(gpu_kmers_valid_p, 0, KMERS_COMBINATIONS * sizeof(qc_kmers_t));
                memset(gpu_kmers_invalid_p, 0, KMERS_COMBINATIONS * sizeof(qc_kmers_t));
            }

            // calculation in CPU
            LOG_DEBUG("qc calculation on CPU thread...\n");

            mean_batch_size += fastq_batch_list_item_p->batch_p->data_size;

            //if (time_flag) { start_timer(t1_cpu); }

            // if trim and filter quality is the same we only calculate it once
            if ((input_p->ltrim_nts == input_p->lfilter_nts) && (input_p->rtrim_nts == input_p->rfilter_nts)) {
                cpu_prepro_trim(fastq_batch_list_item_p->batch_p->data_indices, fastq_batch_list_item_p->batch_p->seq, fastq_batch_list_item_p->batch_p->quality, gpu_result_p, fastq_batch_list_item_p->batch_p->num_reads, input_p->min_quality, input_p->max_quality, input_p->rtrim_nts, input_p->ltrim_nts, input_p->begin_quality_nt, input_p->end_quality_nt);
            } else {
                cpu_prepro_filter_trim(fastq_batch_list_item_p->batch_p->data_indices, fastq_batch_list_item_p->batch_p->seq, fastq_batch_list_item_p->batch_p->quality, gpu_result_p, fastq_batch_list_item_p->batch_p->num_reads, input_p->min_quality, input_p->max_quality, input_p->rtrim_nts, input_p->ltrim_nts, input_p->rfilter_nts, input_p->lfilter_nts, input_p->begin_quality_nt, input_p->end_quality_nt);
            }

            //if (time_flag) { stop_timer(t1_cpu, t2_cpu, cpu_time); }

            // create a new qc_batch object
            qc_batch_p = (qc_batch_t*) malloc(sizeof(qc_batch_t));
            qc_batch_p->id = fastq_batch_list_item_p->id;
            qc_batch_p->source_id = fastq_batch_list_item_p->batch_p->source_id;
            qc_batch_p->nb_reads = fastq_batch_list_item_p->batch_p->num_reads;
            qc_batch_p->gpu_result_p = gpu_result_p;
            qc_batch_p->gpu_kmers_valid_p = gpu_kmers_valid_p;
            qc_batch_p->gpu_kmers_invalid_p = gpu_kmers_invalid_p;
            qc_batch_p->gpu_nt_type_valid_counter_p = NULL;
            qc_batch_p->gpu_nt_type_invalid_counter_p = NULL;
            qc_batch_p->read_p = fastq_batch_list_item_p->batch_p;

            // insert batch into the list
            item_p = list_item_new(qc_batch_p->id, 0, qc_batch_p);
            list_insert_item(item_p, &qc_batch_list);

            // free the current batch item
            fastq_batch_list_item_free(fastq_batch_list_item_p, 0);

            sprintf(log_message, "Thread-GPU:...processing for batch %i done !\n", qc_batch_p->id);
            LOG_DEBUG(log_message);
        } // end if-else

        // ask again for reads server status
        reads_alive = fastq_batch_list_get_producers(batch_list_p);

        // next batch: the first in the list
        fastq_batch_list_item_p = fastq_batch_list_remove(batch_list_p);
    } // end of loop

    pthread_mutex_lock(&gpus_thread_alive_lock);
    gpus_thread_alive--;
    pthread_mutex_unlock(&gpus_thread_alive_lock);

    list_decr_writers(&qc_batch_list);

    LOG_DEBUG("Thread-GPU: END\n");

    // exit pthread
    pthread_exit(0);
}

/* **********************************************
 *    		Results server   		*
 * *********************************************/

void* results_server(void* params_p) {
    LOG_DEBUG("Thread-RESULTS: START\n");
    int source_id = 0;
    results_server_input_t* input_p = (results_server_input_t*) params_p;

    int reads, length, mean_read_quality;
    int base_quality = input_p->base_quality;
    int cpu_num_threads = input_p->cpu_num_threads;
    int cg_flag = input_p->cg_flag;
    int k_cg = input_p->k_cg;
    char* gs_filename = input_p->genomic_signature_input;

    //vectors for read sum length
    long read_sum_length_valid[MAX_NUM_PRODUCERS];
    memset(&read_sum_length_valid, 0, sizeof(read_sum_length_valid));
    long read_sum_length_invalid[MAX_NUM_PRODUCERS];
    memset(&read_sum_length_invalid, 0, sizeof(read_sum_length_invalid));

    //data initialization for chaos game (gs: genomic signature)
    chaos_game_data_t** chaos_game_data_valid_p;
    chaos_game_data_t** chaos_game_data_invalid_p;
    unsigned int** table_gs_p;
    header_gs_t* header_gs_p;
    unsigned int ref_word_count;
    int dim_n = (1 << k_cg); //dim_n = 1 << k_cg (2^k)

    if (cg_flag) {
        if (time_flag) {
            start_timer(t1_chaos_game);
        }

        chaos_game_data_valid_p = (chaos_game_data_t**) calloc(input_p->num_sources, sizeof(chaos_game_data_t*));
        chaos_game_data_invalid_p = (chaos_game_data_t**) calloc(input_p->num_sources, sizeof(chaos_game_data_t*));
        table_gs_p = chaos_game_table_gs_init(table_gs_p, dim_n);
        header_gs_p = (header_gs_t*) calloc(1, sizeof(header_gs_t));

        ref_word_count = chaos_game_load_table_gs_direct(table_gs_p, dim_n, gs_filename);

        for (int n = 0; n < input_p->num_sources; n++) {
            chaos_game_data_valid_p[n] = (chaos_game_data_t*) calloc(1, sizeof(chaos_game_data_t));
            chaos_game_data_invalid_p[n] = (chaos_game_data_t*) calloc(1, sizeof(chaos_game_data_t));
            chaos_game_data_init(chaos_game_data_valid_p[n], table_gs_p, k_cg, dim_n, input_p->base_quality);
            chaos_game_data_init(chaos_game_data_invalid_p[n], table_gs_p, k_cg, dim_n, input_p->base_quality);
        }

        header_gs_init(header_gs_p, gs_filename, k_cg, chaos_game_data_valid_p[0]->memory_size);

        if (time_flag) {
            stop_timer(t1_chaos_game, t2_chaos_game, chaos_game_time);
        }
    }

    //qc_report declarations and memory allocation
    qc_report_t* qc_report_valid = (qc_report_t*) calloc(MAX_NUM_PRODUCERS, sizeof(qc_report_t));
    qc_report_t* qc_report_invalid = (qc_report_t*) calloc(MAX_NUM_PRODUCERS, sizeof(qc_report_t));

    // go through the results batch list, and process it
    // take and remove the first item, and so on...
    qc_batch_t* qc_batch_p = NULL;
    list_item_t* item_p = NULL;
    int gpus_alive;

    // variables for quality per nucleotide calculation
    if ((input_p->qc_step) && (input_p->kmers_on)) {
        char kmer[6];
        for (int k = 0; k < input_p->num_sources; k++) {
            qc_report_valid[k].min_read_length = 10000000;
            qc_report_valid[k].max_read_length = 0;
            qc_report_invalid[k].min_read_length = 10000000;
            qc_report_invalid[k].max_read_length = 0;

            if (input_p->kmers_on) {
                qc_report_valid[k].qc_kmers = (qc_kmers_t*) malloc(KMERS_COMBINATIONS * sizeof(qc_kmers_t));
                qc_report_invalid[k].qc_kmers = (qc_kmers_t*) malloc(KMERS_COMBINATIONS * sizeof(qc_kmers_t));

                for (int j = 0; j < KMERS_COMBINATIONS; j++) {
                    qc_report_valid[k].qc_kmers[j].id = j;
                    memcpy(&qc_report_valid[k].qc_kmers[j].kmer, kmers_string(j, kmer), 5 * sizeof(char));
                    qc_report_valid[k].qc_kmers[j].total_count = 0;
                    memset(qc_report_valid[k].qc_kmers[j].position_count, 0, MAX_LINE_LENGTH * sizeof(int));

                    qc_report_invalid[k].qc_kmers[j].id = j;
                    memcpy(&qc_report_invalid[k].qc_kmers[j].kmer, kmers_string(j, kmer), 5 * sizeof(char));
                    qc_report_invalid[k].qc_kmers[j].total_count = 0;
                    memset(qc_report_invalid[k].qc_kmers[j].position_count, 0, MAX_LINE_LENGTH * sizeof(int));
                }
            }
        }
    }

    pthread_mutex_lock(&gpus_thread_alive_lock);
    gpus_alive = gpus_thread_alive;
    pthread_mutex_unlock(&gpus_thread_alive_lock);

    while (gpus_alive > 0 || qc_batch_p != NULL) {
        if (qc_batch_p == NULL) {
            // Delay for a bit
            sched_yield();
            usleep(1000);
        } else {
            number_of_batchs++;

            if (time_flag) {
                start_timer(t1_result);
            }

            source_id = qc_batch_p->source_id;

            char log_message[100];
            sprintf(log_message, "Thread-RESULTS: processing for batch %i....\n", qc_batch_p->id);
            LOG_DEBUG(log_message);

            // result processing batch per batch
            reads = qc_batch_p->nb_reads;
            mean_reads_per_batch += reads;

            qc_batch_p->gpu_nt_type_valid_counter_p = (int*) malloc(COUNTERS_SIZE * MAX_LINE_LENGTH * sizeof(int));
            memset(qc_batch_p->gpu_nt_type_valid_counter_p, 0, COUNTERS_SIZE * MAX_LINE_LENGTH * sizeof(int));
            qc_batch_p->gpu_nt_type_invalid_counter_p = (int*) malloc(COUNTERS_SIZE * MAX_LINE_LENGTH * sizeof(int));
            memset(qc_batch_p->gpu_nt_type_invalid_counter_p, 0, COUNTERS_SIZE * MAX_LINE_LENGTH * sizeof(int));

            status_batch_t* status_batch_p = (status_batch_t*) malloc(sizeof(status_batch_t));
            status_batch_p->read_status = (char*) malloc(reads * sizeof(char));
            status_batch_p->id = qc_batch_p->id;
            status_batch_p->source_id = qc_batch_p->source_id;
            status_batch_p->num_reads = reads;
            status_batch_p->read_p = qc_batch_p->read_p;
            status_batch_p->prev_p = NULL;
            status_batch_p->next_p = NULL;

            // calculation of read status based on: read length, read quality, no-determined (N) positions and last nt mean quality
            if (time_flag) {
                stop_timer(t1_result, t2_result, result_time);
            }
            if (time_flag) {
                start_timer(t1_cpu);
            }
            
            for (int i = 0 ; i < reads ; i++) {
                length = qc_batch_p->gpu_result_p[i].counters[TOTAL];

                status_batch_p->read_status[i] = UNKNOWN_STATUS_READ;
                if (input_p->prepro_step) {
                    preprocessing_read(i, input_p, qc_batch_p, status_batch_p);
                }

                if (input_p->filter_step) {
                    filtering_read(i, length, input_p, qc_batch_p, status_batch_p);
                }
                //if (time_flag) { start_timer(t1_cpu); }
                if (input_p->qc_step) {
                    if (status_batch_p->read_status[i] == INVALID_READ) {
                        qc_report_invalid[source_id].nb_reads++;
                        qc_per_read(i, length, base_quality, &mean_read_quality, read_sum_length_invalid, input_p, qc_batch_p, status_batch_p, qc_report_invalid, cpu_num_threads);
                    } else {
                        qc_report_valid[source_id].nb_reads++;
                        qc_per_read(i, length, base_quality, &mean_read_quality, read_sum_length_valid, input_p, qc_batch_p, status_batch_p, qc_report_valid, cpu_num_threads);
                    }
                }
                //if (time_flag) { stop_timer(t1_cpu, t2_cpu, cpu_time); }
            }

            if (time_flag) {
                stop_timer(t1_cpu, t2_cpu, cpu_time);
            }
            if (time_flag) {
                start_timer(t1_kmers);
            }

            if (input_p->kmers_on) {
                pthread_mutex_t kmers_lock = PTHREAD_MUTEX_INITIALIZER;
                int num_thread;
                int** kmers_valid_thread_counter = (int**) calloc(cpu_num_threads, sizeof(int*));
                int** kmers_invalid_thread_counter = (int**) calloc(cpu_num_threads, sizeof(int*));

                for (int n = 0; n < cpu_num_threads; n++) {
                    kmers_valid_thread_counter[n] = (int*) calloc(KMERS_COMBINATIONS, sizeof(int));
                    kmers_invalid_thread_counter[n] = (int*) calloc(KMERS_COMBINATIONS, sizeof(int));
                }

                #pragma omp parallel for num_threads(cpu_num_threads) shared(kmers_lock, kmers_valid_thread_counter, kmers_invalid_thread_counter) schedule(dynamic, 1000)
                for (int i = 0 ; i < reads ; i++) {
                    num_thread = omp_get_thread_num();
                    length = qc_batch_p->gpu_result_p[i].counters[TOTAL];

                    if (length > 0) {
                        kmers_thread_per_read(i, length, input_p, qc_batch_p, kmers_valid_thread_counter[num_thread], kmers_invalid_thread_counter[num_thread], status_batch_p, cpu_num_threads);
                    }
                }
                #pragma omp end parallel

                //reduction/accumulation of the partial count kmers vectors
                for (int n = 0; n < cpu_num_threads; n++) {
                    for (int p = 0; p < KMERS_COMBINATIONS; p++) {
                        qc_batch_p->gpu_kmers_valid_p[p].total_count += kmers_valid_thread_counter[n][p];
                        qc_batch_p->gpu_kmers_invalid_p[p].total_count += kmers_invalid_thread_counter[n][p];
                    }
                }

                //free counters memory
                for (int n = 0; n < cpu_num_threads; n++) {
                    free(kmers_valid_thread_counter[n]);
                    free(kmers_invalid_thread_counter[n]);
                }
                free(kmers_valid_thread_counter);
                free(kmers_invalid_thread_counter);
            }

            if (time_flag) {
                stop_timer(t1_kmers, t2_kmers, kmers_time);
            }
            printf("qc_report_invalid[0].nb_reads: %lu, qc_report_valid[0].nb_reads: %lu\n", qc_report_invalid[0].nb_reads, qc_report_valid[0].nb_reads);

            if (time_flag) {
                start_timer(t1_result);
            }

            if (input_p->qc_step) {
                if (time_flag) {
                    stop_timer(t1_result, t2_result, result_time);
                }

                qc_per_batch(source_id, input_p, qc_batch_p, qc_report_invalid, 0);
                qc_per_batch(source_id, input_p, qc_batch_p, qc_report_valid, 1);

                if (time_flag) {
                    start_timer(t1_result);
                }
            }

            //fill chaos game table sequence table batch by batch
            if ((input_p->cg_flag) && (qc_batch_p->id < 2)) {
                if (time_flag) {
                    start_timer(t1_chaos_game);
                }

                chaos_game_fill_tables(chaos_game_data_valid_p[source_id], status_batch_p, ALL_READS);

                if (time_flag) {
                    stop_timer(t1_chaos_game, t2_chaos_game, chaos_game_time);
                }
            }

            sprintf(log_message, "Thread-RESULTS: ...processing for batch %i done!\n", qc_batch_p->id);
            LOG_DEBUG(log_message);

            // if preprocessing of filtering is enable reads are maintained in memory and status batch is inserted in the list
            // otherwise we free all resources since writer thread is not going to be run
            if ((input_p->prepro_step) || (input_p->filter_step)) {
                status_batch_list_insert(status_batch_p, &status_batch_list);
                if (qc_batch_p->gpu_result_p != NULL) {
                    free(qc_batch_p->gpu_result_p);
                }
                qc_batch_free(qc_batch_p, 0);
            } else {
                int sid = qc_batch_p->source_id;
                qc_batch_free(qc_batch_p, 1);
                status_batch_free(status_batch_p, 1);
            }

            if (time_flag) {
                stop_timer(t1_result, t2_result, result_time);
            }
        } //end of if-else

        // getting gpus thread status
        pthread_mutex_lock(&gpus_thread_alive_lock);
        gpus_alive = gpus_thread_alive;
        pthread_mutex_unlock(&gpus_thread_alive_lock);

        // next batch
        list_item_free(item_p);
        item_p = list_remove_item(&qc_batch_list);

        if (item_p != NULL) {
            qc_batch_p = (qc_batch_t*) item_p->data_p;
        } else {
            break;
        }
    } // end of batch loop

    // after chaos game calculations its structures must be freed
    // qc_report is also filled with the necesary data for the report
    if ((input_p->qc_step) && (cg_flag)) {
        if (time_flag) {
            start_timer(t1_chaos_game);
        }

        chaos_game_load_table_gs(chaos_game_data_valid_p[0], gs_filename); //load genomic signature for reference sequence only once

        for (int n = 0; n < input_p->num_sources; n++) {
            //chaos_game_print_table(chaos_game_data_valid_p[n]->table_q, chaos_game_data_valid_p[n]->dim_n);
            chaos_game_data_valid_p[n]->table_gs = chaos_game_data_valid_p[0]->table_gs;
            chaos_game_data_valid_p[n]->ref_word_count = chaos_game_data_valid_p[0]->ref_word_count;
	    
            chaos_game_calculate_table_dif(chaos_game_data_valid_p[n]);  //calculate diff table between fastq sequences and the reference

            if (chaos_game_validate_table_dif(chaos_game_data_valid_p[n])) { //if diff table is not valid then generate files and images
                chaos_game_write_table_images(chaos_game_data_valid_p[n], input_p->source_p[n].filename, input_p->report_directory);
            }

            qc_report_get_values_from_chaos_game(chaos_game_data_valid_p[n], &qc_report_valid[n]);
        }

        for (int n = 0; n < input_p->num_sources; n++) {
            chaos_game_data_free(chaos_game_data_valid_p[n]);
            chaos_game_data_free(chaos_game_data_invalid_p[n]);
        }

        free(chaos_game_data_valid_p);
        free(chaos_game_data_invalid_p);
        chaos_game_table_gs_free(table_gs_p, dim_n);
        free(header_gs_p);
        if (time_flag) {
            stop_timer(t1_chaos_game, t2_chaos_game, chaos_game_time);
        }
    }

    // accumulate final results and generate report
    if (input_p->qc_step) {
        qc_results(input_p, read_sum_length_valid, qc_report_valid, 1);

        if (input_p->filter_step) {
            qc_results(input_p, read_sum_length_invalid, qc_report_invalid, 0);
        }
    }

    // qc_report structures must be freed
    if (input_p->kmers_on) {
        for (int k = 0; k < input_p->num_sources; k++) {
            free(qc_report_valid[k].qc_kmers);
            free(qc_report_invalid[k].qc_kmers);
        }
    }
    free(qc_report_valid);
    free(qc_report_invalid);

    pthread_mutex_lock(&results_thread_alive_lock);
    results_thread_alive--;
    pthread_mutex_unlock(&results_thread_alive_lock);

    LOG_DEBUG("Thread-RESULTS: END\n");

    // exiting pthread
    pthread_exit(0);
}

/* ******************************************************
 *    		Writer single-end server   		*
 * *****************************************************/

void* writer_single_end_server(void* params_p) {
    LOG_DEBUG("Thread-WRITER: START\n");

    int source_id = 0;
    writer_server_input_t* input_p = (writer_server_input_t*) params_p;

    int reads;
    char input_shortname[MAX_FULL_PATH_LENGTH];
    char output_filename[MAX_FULL_PATH_LENGTH];

    //file descriptors for preprocessed reads
    get_filename_from_path(input_p->source_p[source_id].filename, input_shortname);
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, VALID_READS_FILE_SUFFIX);
    FILE* fd_valid = fopen(output_filename, "w");
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, INVALID_READS_FILE_SUFFIX);
    FILE* fd_invalid = fopen(output_filename, "w");

    // go through the status_batch list, process it and write the valid and invalid reads back
    status_batch_t* status_batch_p = NULL;
    int results_alive;

    // getting results thread status
    pthread_mutex_lock(&results_thread_alive_lock);
    results_alive = results_thread_alive;
    pthread_mutex_unlock(&results_thread_alive_lock);

    // get the first element in the list
    status_batch_p = status_batch_list_remove(&status_batch_list);

    while (results_alive > 0 || status_batch_p != NULL) {
        if (status_batch_p == NULL) {
            // Delay for a bit
            sched_yield();
            usleep(1000);
        } else {
            if (time_flag) {
                start_timer(t1_write);
            }

            source_id = status_batch_p->source_id;

            char log_message[100];
            sprintf(log_message, "Thread-WRITER: processing for status batch %i....\n", status_batch_p->id);
            LOG_DEBUG(log_message);

            // read processing batch per batch
            reads = status_batch_p->num_reads;

            // calculation of read status based on: read length, read quality, no-determined (N) positions and last nt mean quality
            for (int i = 0 ; i < reads ; i++) {
                if (status_batch_p->read_status[i] == (char) VALID_READ) {
                    fprintf_read(fd_valid, status_batch_p->read_p, i);
                } else if (status_batch_p->read_status[i] == (char) INVALID_READ) {
                    fprintf_read(fd_invalid, status_batch_p->read_p, i);
                } else if (status_batch_p->read_status[i] == (char) RTRIM_READ) {
                    fprintf_rtrim_read(fd_valid, status_batch_p->read_p, i, input_p->rtrim_nts);
                } else if (status_batch_p->read_status[i] == (char) TRIM_READ) {
                    fprintf_trim_read(fd_valid, status_batch_p->read_p, i, input_p->rtrim_nts, input_p->ltrim_nts);
                } else if (status_batch_p->read_status[i] == (char) LTRIM_READ) {
                    fprintf_ltrim_read(fd_valid, status_batch_p->read_p, i, input_p->ltrim_nts);
                } else {
                    LOG_ERROR("Invalid read status, not VALID_READ, INVALID_READ or CUT_READ");
                }
            }

            // free ALL memory
            status_batch_free(status_batch_p, 1);

            if (time_flag) {
                stop_timer(t1_write, t2_write, write_time);
            }

            sprintf(log_message, "Thread-WRITER: ....processing for status batch %i done !\n", status_batch_p->id);
            LOG_DEBUG(log_message);
        } //end of if-else

        // getting gpus thread status
        pthread_mutex_lock(&results_thread_alive_lock);
        results_alive = results_thread_alive;
        pthread_mutex_unlock(&results_thread_alive_lock);

        // next batch
        status_batch_p = status_batch_list_remove(&status_batch_list);
    } // end of batch loop

    fclose(fd_valid);
    fclose(fd_invalid);

    LOG_DEBUG("Thread-WRITER: END\n");

    // exiting pthread
    pthread_exit(0);
}

/* ******************************************************
 *    		Writer paired-end server   		*
 * *****************************************************/

void* writer_paired_end_server(void* params_p) {
    LOG_DEBUG("Thread-WRITER: START\n");
    int source_id = 0;
    writer_server_input_t* input_p = (writer_server_input_t*) params_p;

    int reads_1, reads_2;
    int processed_reads_1 = 0;
    int processed_reads_2 = 0;
    char input_shortname[MAX_FULL_PATH_LENGTH];
    char output_filename[MAX_FULL_PATH_LENGTH];

    //file descriptors for reads preprocessed
    get_filename_from_path(input_p->source_p[0].filename, input_shortname);
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, VALID_READS_FILE_SUFFIX);
    FILE* fd_valid_1 = fopen(output_filename, "w");
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, INVALID_READS_FILE_SUFFIX);
    FILE* fd_invalid_1 = fopen(output_filename, "w");
    get_filename_from_path(input_p->source_p[1].filename, input_shortname);
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, VALID_READS_FILE_SUFFIX);
    FILE* fd_valid_2 = fopen(output_filename, "w");
    sprintf(output_filename, "%s/%s%s", input_p->output_directory, input_shortname, INVALID_READS_FILE_SUFFIX);
    FILE* fd_invalid_2 = fopen(output_filename, "w");

    // go through the status_batch list, process it and write the valid and invalid reads back...
    status_batch_t* status_batch_1_p = NULL;
    status_batch_t* status_batch_2_p = NULL;
    int results_alive;

    int batch1_count = 0;
    int batch2_count = 0;

    // getting results thread status
    pthread_mutex_lock(&results_thread_alive_lock);
    results_alive = results_thread_alive;
    pthread_mutex_unlock(&results_thread_alive_lock);

    // get the first element in the list
    status_batch_1_p = status_batch_list_get_next_by_id(&status_batch_list, input_p->source_p[0].id, batch1_count);
    status_batch_2_p = status_batch_list_get_next_by_id(&status_batch_list, input_p->source_p[1].id, batch2_count);

    while (results_alive > 0 || status_batch_1_p != NULL || status_batch_2_p != NULL) {
        if (status_batch_1_p == NULL || status_batch_2_p == NULL) {
            // Delay for a bit
            sched_yield();
            usleep(100);
        } else {
            if (time_flag) {
                start_timer(t1_write);
            }

            char log_message[100];
            sprintf(log_message, "Thread-WRITER: processing for status batch 1: %i and status batch 2: %i....\n", status_batch_1_p->id, status_batch_2_p->id);
            LOG_DEBUG(log_message);

            // read processing batch per batch
            reads_1 = status_batch_1_p->num_reads;
            reads_2 = status_batch_2_p->num_reads;

            // calculation of read status based on: read length, read quality, no-determined (N) positions and last nt mean quality
            while ((processed_reads_1 < reads_1) && (processed_reads_2 < reads_2)) {
                char* find = strchr(&(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]), ' ');
                if (find == NULL) {
                    find = strchr(&(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]), '\t');
                    if (find == NULL) {
                        find = strchr(&(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]), '\0');
                    }
                }

                int length = find - &(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]);

                if (strncmp(&(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]), &(status_batch_2_p->read_p->header[status_batch_2_p->read_p->header_indices[processed_reads_2]]), length) == 0) {
                    if ((status_batch_1_p->read_status[processed_reads_1] == (char) VALID_READ) && (status_batch_2_p->read_status[processed_reads_2] == (char) VALID_READ)) {
                        fprintf_read(fd_valid_1, status_batch_1_p->read_p, processed_reads_1);
                        fprintf_read(fd_valid_2, status_batch_2_p->read_p, processed_reads_2);
                    } else if ((status_batch_1_p->read_status[processed_reads_1] == (char) INVALID_READ) || (status_batch_2_p->read_status[processed_reads_2] == (char) INVALID_READ)) {
                        fprintf_read(fd_invalid_1, status_batch_1_p->read_p, processed_reads_1);
                        fprintf_read(fd_invalid_2, status_batch_2_p->read_p, processed_reads_2);
                    } else {
                        if (status_batch_1_p->read_status[processed_reads_1] == (char) RTRIM_READ) {
                            fprintf_rtrim_read(fd_valid_1, status_batch_1_p->read_p, processed_reads_1, input_p->rtrim_nts);
                        } else if (status_batch_1_p->read_status[processed_reads_1] == (char) TRIM_READ) {
                            fprintf_trim_read(fd_valid_1, status_batch_1_p->read_p, processed_reads_1, input_p->rtrim_nts, input_p->ltrim_nts);
                        } else if (status_batch_1_p->read_status[processed_reads_1] == (char) LTRIM_READ) {
                            fprintf_ltrim_read(fd_valid_1, status_batch_1_p->read_p, processed_reads_1, input_p->ltrim_nts);
                        } else {
                            fprintf_read(fd_valid_1, status_batch_1_p->read_p, processed_reads_1);
                        }

                        if (status_batch_2_p->read_status[processed_reads_2] == (char) RTRIM_READ) {
                            fprintf_rtrim_read(fd_valid_2, status_batch_2_p->read_p, processed_reads_2, input_p->rtrim_nts);
                        } else if (status_batch_2_p->read_status[processed_reads_2] == (char) TRIM_READ) {
                            fprintf_trim_read(fd_valid_2, status_batch_2_p->read_p, processed_reads_2, input_p->rtrim_nts, input_p->ltrim_nts);
                        } else if (status_batch_2_p->read_status[processed_reads_2] == (char) LTRIM_READ) {
                            fprintf_ltrim_read(fd_valid_2, status_batch_2_p->read_p, processed_reads_2, input_p->ltrim_nts);
                        } else {
                            fprintf_read(fd_valid_2, status_batch_2_p->read_p, processed_reads_2);
                        }
                    }
                } else {
                    fprintf_read(fd_invalid_1, status_batch_1_p->read_p, processed_reads_1);
                    fprintf_read(fd_invalid_2, status_batch_2_p->read_p, processed_reads_2);

                    sprintf(log_message, "Distinct headers!!! %s != %s\n", &(status_batch_1_p->read_p->header[status_batch_1_p->read_p->header_indices[processed_reads_1]]), &(status_batch_2_p->read_p->header[status_batch_2_p->read_p->header_indices[processed_reads_2]]));
                    LOG_DEBUG(log_message);
                }

                processed_reads_1++;
                processed_reads_2++;
            }

            sprintf(log_message, "Thread-WRITER: ....processing for status batch 1: %i and status batch 2: %i done !\n", status_batch_1_p->id, status_batch_2_p->id);
            LOG_DEBUG(log_message);

            // free ALL memory
            if (reads_1 == processed_reads_1) {
                status_batch_free(status_batch_1_p, 1);
                status_batch_1_p = NULL;
                processed_reads_1 = 0;
                batch1_count++;
            }

            if (reads_2 == processed_reads_2) {
                status_batch_free(status_batch_2_p, 1);
                status_batch_2_p = NULL;
                processed_reads_2 = 0;
                batch2_count++;
            }

            if (time_flag) {
                stop_timer(t1_write, t2_write, write_time);
            }
        } //end of if-else

        // next batch...
        if (status_batch_1_p == NULL) {
            status_batch_1_p = status_batch_list_get_next_by_id(&status_batch_list, input_p->source_p[0].id, batch1_count);
        }
        if (status_batch_2_p == NULL) {
            status_batch_2_p = status_batch_list_get_next_by_id(&status_batch_list, input_p->source_p[1].id, batch2_count);
        }

        // getting gpus thread status
        pthread_mutex_lock(&results_thread_alive_lock);
        results_alive = results_thread_alive;
        pthread_mutex_unlock(&results_thread_alive_lock);

    } // end of batch loop

    fclose(fd_valid_1);
    fclose(fd_invalid_1);
    fclose(fd_valid_2);
    fclose(fd_invalid_2);
    LOG_DEBUG("Thread-WRITER: END\n");

    // exiting pthread
    pthread_exit(0);
}

/* **************************************************************************************
 *    		Single-end processing main function implementation (public)  		*
 * *************************************************************************************/

void kernel_prepro_fastq_single_end(size_t batch_size, int max_fastq_batch_list_length, int num_blocks, int num_threads, int cpu_num_threads, int cpu_qc_calc_num_threads, char* filename, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input) {
    int num_gpu_devices = 0;

    gpus_thread_alive = num_gpu_devices + cpu_qc_calc_num_threads;

    int num_total_threads = num_blocks * num_threads;

    // some local variables
    void* r;

    // multi-threads, one gpu_server_thread by gpu device
    pthread_t** qc_calc_server_thread_p = (pthread_t**) calloc((num_gpu_devices + cpu_qc_calc_num_threads), sizeof(pthread_t*));
    qc_calc_server_input_t** qc_calc_server_input_p = (qc_calc_server_input_t**) calloc((num_gpu_devices + cpu_qc_calc_num_threads), sizeof(qc_calc_server_input_t*));

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        qc_calc_server_thread_p[i] = (pthread_t*) calloc(1, sizeof(pthread_t));
        qc_calc_server_input_p[i] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
    }

    pthread_t results_server_thread;
    pthread_t writer_server_thread;

    //only one fastq_batch_list
    fastq_batch_list_t batch_list;
    fastq_batch_list_init(&batch_list, 0);
    status_batch_list_init(&status_batch_list);

    if (num_gpu_devices > 0) {
        list_init("qc_batch_list", (num_gpu_devices + cpu_qc_calc_num_threads), max_fastq_batch_list_length, &qc_batch_list);
    } else {
        list_init("qc_batch_list", cpu_qc_calc_num_threads, max_fastq_batch_list_length, &qc_batch_list);
    }

    // calling thread to serve reads from file, but first, prepare input parameter
    int num_sources = 1;
    source_t source_p;
    source_p.id = 0;
    strcpy(source_p.filename, filename);

    fastq_batch_reader_t* batch_reader_p = fastq_batch_reader_new(filename, source_p.id, &batch_list, batch_size, &qc_batch_list, max_fastq_batch_list_length);

    fastq_batch_reader_start(batch_reader_p);

    // calling thread to serve GPUs, and again, prepare input parameter
    if (num_gpu_devices > 0) { // GPU implementacion
        for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
            qc_calc_server_input_p[i]->num_gpu_devices = num_gpu_devices;
            qc_calc_server_input_p[i]->cpu_num_threads = 0;
            qc_calc_server_input_p[i]->gpu_device_id[0] = i;
            qc_calc_server_input_p[i]->num_blocks = num_blocks;
            qc_calc_server_input_p[i]->num_threads = num_threads;
            qc_calc_server_input_p[i]->min_quality = min_quality;
            qc_calc_server_input_p[i]->max_quality = max_quality;
            qc_calc_server_input_p[i]->begin_quality_nt = begin_quality_nt;
            qc_calc_server_input_p[i]->end_quality_nt = end_quality_nt;
            qc_calc_server_input_p[i]->rtrim_nts = rtrim_nts;
            qc_calc_server_input_p[i]->ltrim_nts = ltrim_nts;
            qc_calc_server_input_p[i]->rfilter_nts = rfilter_nts;
            qc_calc_server_input_p[i]->lfilter_nts = lfilter_nts;
            qc_calc_server_input_p[i]->kmers_on = kmers_on;
            qc_calc_server_input_p[i]->batch_list_p = &batch_list;
            pthread_create(qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
        }
    } else {
        for (int i = 0; i < cpu_qc_calc_num_threads; i++) {
            qc_calc_server_input_p[i]->num_gpu_devices = 0;
            qc_calc_server_input_p[i]->cpu_num_threads = cpu_num_threads;
            qc_calc_server_input_p[i]->gpu_device_id[0] = 0;
            qc_calc_server_input_p[i]->num_blocks = num_blocks;
            qc_calc_server_input_p[i]->num_threads = num_threads;
            qc_calc_server_input_p[i]->min_quality = min_quality;
            qc_calc_server_input_p[i]->max_quality = max_quality;
            qc_calc_server_input_p[i]->begin_quality_nt = begin_quality_nt;
            qc_calc_server_input_p[i]->end_quality_nt = end_quality_nt;
            qc_calc_server_input_p[i]->rtrim_nts = rtrim_nts;
            qc_calc_server_input_p[i]->ltrim_nts = ltrim_nts;
            qc_calc_server_input_p[i]->rfilter_nts = rfilter_nts;
            qc_calc_server_input_p[i]->lfilter_nts = lfilter_nts;
            qc_calc_server_input_p[i]->kmers_on = kmers_on;
            qc_calc_server_input_p[i]->batch_list_p = &batch_list;
            pthread_create(qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
        }
    }

    // calling thread to process results from QC-CALC (GPU or CPU)
    results_server_input_t results_server_input;
    results_server_input.num_blocks = num_blocks;
    results_server_input.num_threads = num_threads;
    results_server_input.num_sources = 1;
    results_server_input.source_p = &source_p;
    results_server_input.report_directory = output_directory;
    results_server_input.min_quality = min_quality;
    results_server_input.max_quality = max_quality;
    results_server_input.base_quality = base_quality;
    results_server_input.max_nts_mismatch =  max_nts_mismatch;
    results_server_input.max_n_per_read = max_n_per_read;
    results_server_input.min_read_length = min_read_length;
    results_server_input.max_read_length = max_read_length;
    results_server_input.rtrim_nts = rtrim_nts;
    results_server_input.ltrim_nts = ltrim_nts;
    results_server_input.rfilter_nts = rfilter_nts;
    results_server_input.lfilter_nts = lfilter_nts;
    results_server_input.prepro_step = prepro_step;
    results_server_input.filter_step = filter_step;
    results_server_input.qc_step = qc_step;
    results_server_input.kmers_on = kmers_on;
    results_server_input.cpu_num_threads =  cpu_num_threads;
    results_server_input.cg_flag = cg_flag;
    results_server_input.k_cg = k_cg;
    results_server_input.genomic_signature_input = genomic_signature_input;
    pthread_create(&results_server_thread, NULL, results_server, (void*) &results_server_input);

    // calling thread to write the preprocessed reads back
    writer_server_input_t writer_server_input;
    writer_server_input.rtrim_nts = rtrim_nts;
    writer_server_input.ltrim_nts = ltrim_nts;
    writer_server_input.num_sources = num_sources;
    writer_server_input.source_p = &source_p;
    writer_server_input.output_directory = output_directory;
    pthread_create(&writer_server_thread, NULL, writer_single_end_server, (void*) &writer_server_input);

    // wait for all terminating
    fastq_batch_reader_join(batch_reader_p);

    pthread_join(results_server_thread, &r);

    if ((prepro_step) || (filter_step)) {
        pthread_join(writer_server_thread, &r);
    }

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        pthread_join(*qc_calc_server_thread_p[i], &r);
        free(qc_calc_server_thread_p[i]);
    }

    // free thread stuff and parameters
    fastq_batch_reader_free(batch_reader_p);

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        free(qc_calc_server_input_p[i]);
    }
    free(qc_calc_server_input_p);
}

/* **************************************************************************************
 *    		Paired-end processing main function implementations (public)  		*
 * *************************************************************************************/

void kernel_prepro_fastq_paired_end(size_t batch_size, int max_fastq_batch_list_length, int num_blocks, int num_threads, int cpu_num_threads, int cpu_qc_calc_num_threads, char* filename1, char* filename2, char* output_directory, int min_quality, int max_quality, int base_quality, int begin_quality_nt, int end_quality_nt, int max_nts_mismatch, int max_n_per_read, int min_read_length, int max_read_length, int rtrim_nts, int ltrim_nts, int rfilter_nts, int lfilter_nts, int prepro_step, int filter_step, int qc_step, int kmers_on, int cg_flag, int k_cg, char* genomic_signature_input) {
    int num_gpu_devices = 0;

    gpus_thread_alive = num_gpu_devices + cpu_qc_calc_num_threads;

    int num_total_threads = num_blocks * num_threads;

    // some local variables
    void* r;

    // multi-threads, one gpu_server_thread by gpu device
    pthread_t** qc_calc_server_thread_p = (pthread_t**) calloc((num_gpu_devices + cpu_qc_calc_num_threads), sizeof(pthread_t*));
    qc_calc_server_input_t** qc_calc_server_input_p = (qc_calc_server_input_t**) calloc((num_gpu_devices + cpu_qc_calc_num_threads), sizeof(qc_calc_server_input_t*));

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        qc_calc_server_thread_p[i] = (pthread_t*) calloc(1, sizeof(pthread_t));
        qc_calc_server_input_p[i] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
    }

    pthread_t results_server_thread;
    pthread_t writer_server_thread;

    //only one fastq_batch_list
    fastq_batch_list_t batch_list;
    fastq_batch_list_init(&batch_list, 0);
    status_batch_list_init(&status_batch_list);

    if (num_gpu_devices > 0) {
        list_init("qc_batch_list", (num_gpu_devices + cpu_qc_calc_num_threads), max_fastq_batch_list_length, &qc_batch_list);
    } else {
        list_init("qc_batch_list", cpu_qc_calc_num_threads, max_fastq_batch_list_length, &qc_batch_list);
    }

    int num_sources = 2;
    source_t* source_p = (source_t*) calloc(num_sources, sizeof(source_t));
    source_p[0].id = 0;
    strcpy(source_p[0].filename, filename1);
    source_p[1].id = 1;
    strcpy(source_p[1].filename, filename2);

    // calling thread to serve reads from file, but first, prepare input parameter
    fastq_batch_reader_t* batch_reader_1_p = fastq_batch_reader_new(source_p[0].filename, source_p[0].id, &batch_list, batch_size, &qc_batch_list, max_fastq_batch_list_length);
    fastq_batch_reader_t* batch_reader_2_p = fastq_batch_reader_new(source_p[1].filename, source_p[1].id, &batch_list, batch_size, &qc_batch_list, max_fastq_batch_list_length);

    if (time_flag) {
        start_timer(t1_read);
    }
    fastq_batch_reader_start(batch_reader_1_p);
    fastq_batch_reader_start(batch_reader_2_p);

    // calling thread to serve GPUs or CPUs, and again, prepare input parameters
    // one thread by gpu
    if (num_gpu_devices > 0) { // GPU implementacion
        for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
            qc_calc_server_input_p[i]->num_gpu_devices = num_gpu_devices;
            qc_calc_server_input_p[i]->cpu_num_threads = 0;
            qc_calc_server_input_p[i]->gpu_device_id[0] = i;
            qc_calc_server_input_p[i]->num_blocks = num_blocks;
            qc_calc_server_input_p[i]->num_threads = num_threads;
            qc_calc_server_input_p[i]->min_quality = min_quality;
            qc_calc_server_input_p[i]->max_quality = max_quality;
            qc_calc_server_input_p[i]->begin_quality_nt = begin_quality_nt;
            qc_calc_server_input_p[i]->end_quality_nt = end_quality_nt;
            qc_calc_server_input_p[i]->rtrim_nts = rtrim_nts;
            qc_calc_server_input_p[i]->ltrim_nts = ltrim_nts;
            qc_calc_server_input_p[i]->rfilter_nts = rfilter_nts;
            qc_calc_server_input_p[i]->lfilter_nts = lfilter_nts;
            qc_calc_server_input_p[i]->kmers_on = kmers_on;
            qc_calc_server_input_p[i]->batch_list_p = &batch_list;
            pthread_create(qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
        }
    } else {
        for (int i = 0; i < cpu_qc_calc_num_threads; i++) {
            qc_calc_server_input_p[i]->num_gpu_devices = 0;
            qc_calc_server_input_p[i]->cpu_num_threads = cpu_num_threads;
            qc_calc_server_input_p[i]->gpu_device_id[0] = 0;
            qc_calc_server_input_p[i]->num_blocks = num_blocks;
            qc_calc_server_input_p[i]->num_threads = num_threads;
            qc_calc_server_input_p[i]->min_quality = min_quality;
            qc_calc_server_input_p[i]->max_quality = max_quality;
            qc_calc_server_input_p[i]->begin_quality_nt = begin_quality_nt;
            qc_calc_server_input_p[i]->end_quality_nt = end_quality_nt;
            qc_calc_server_input_p[i]->rtrim_nts = rtrim_nts;
            qc_calc_server_input_p[i]->ltrim_nts = ltrim_nts;
            qc_calc_server_input_p[i]->rfilter_nts = rfilter_nts;
            qc_calc_server_input_p[i]->lfilter_nts = lfilter_nts;
            qc_calc_server_input_p[i]->kmers_on = kmers_on;
            qc_calc_server_input_p[i]->batch_list_p = &batch_list;
            pthread_create(qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
        }
    }

    // calling thread to process results from QC-CALC (GPU or CPU)
    results_server_input_t results_server_input;
    results_server_input.num_blocks = num_blocks;
    results_server_input.num_threads = num_threads;
    results_server_input.num_sources = num_sources;
    results_server_input.source_p = source_p;
    results_server_input.report_directory = output_directory;
    results_server_input.min_quality = min_quality;
    results_server_input.max_quality = max_quality;
    results_server_input.base_quality = base_quality;
    results_server_input.max_nts_mismatch =  max_nts_mismatch;
    results_server_input.max_n_per_read = max_n_per_read;
    results_server_input.min_read_length = min_read_length;
    results_server_input.max_read_length = max_read_length;
    results_server_input.rtrim_nts = rtrim_nts;
    results_server_input.ltrim_nts = ltrim_nts;
    results_server_input.rfilter_nts = rfilter_nts;
    results_server_input.lfilter_nts = lfilter_nts;
    results_server_input.prepro_step = prepro_step;
    results_server_input.filter_step = filter_step;
    results_server_input.qc_step = qc_step;
    results_server_input.kmers_on = kmers_on;
    results_server_input.cpu_num_threads =  cpu_num_threads;
    results_server_input.cg_flag = cg_flag;
    results_server_input.k_cg = k_cg;
    results_server_input.genomic_signature_input = genomic_signature_input;
    pthread_create(&results_server_thread, NULL, results_server, (void*) &results_server_input);

    // calling thread to write the preprocessed reads back
    writer_server_input_t writer_server_input;
    writer_server_input.rtrim_nts = rtrim_nts;
    writer_server_input.ltrim_nts = ltrim_nts;
    writer_server_input.num_sources = num_sources;
    writer_server_input.source_p = source_p;
    writer_server_input.output_directory = output_directory;

    if ((prepro_step) || (filter_step)) {
        pthread_create(&writer_server_thread, NULL, writer_paired_end_server, (void*) &writer_server_input);
    }

    // wait for all terminating
    // fastq_batch_reader_join returns the number of reads of each file, it must be equal
    unsigned int num_reads_1, num_reads_2;
    num_reads_1 = fastq_batch_reader_join(batch_reader_1_p);
    num_reads_2 = fastq_batch_reader_join(batch_reader_2_p);

    LOG_IF(LOG_FATAL_LEVEL, num_reads_1 != num_reads_2, "Mismatch in number of read of paired end files. Aborting execution.");

    if (time_flag) {
        stop_timer(t1_read, t2_read, read_time);
    }

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        pthread_join(*qc_calc_server_thread_p[i], &r);
        free(qc_calc_server_thread_p[i]);
    }
    free(qc_calc_server_thread_p);

    pthread_join(results_server_thread, &r);

    if ((prepro_step) || (filter_step)) {
        pthread_join(writer_server_thread, &r);
    }

    // free thread stuff and parameters
    free(source_p);
    fastq_batch_reader_free(batch_reader_1_p);
    fastq_batch_reader_free(batch_reader_2_p);

    for (int i = 0; i < (num_gpu_devices + cpu_qc_calc_num_threads); i++) {
        free(qc_calc_server_input_p[i]);
    }
    free(qc_calc_server_input_p);
}

#endif /* PREPRO_CU */
