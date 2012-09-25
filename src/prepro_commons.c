/*
 * Copyright (c) 2012 Victor Requena (BULL)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
 *
 * This file is part of hpg-fastq.
 *
 * hpg-fastq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hpg-fastq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-fastq. If not, see <http://www.gnu.org/licenses/>.
 */


#include "prepro_commons.h"

/* **************************************************************
 *    		Private functions implementations  		*
 * **************************************************************/

void preprocessing_read(int i, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p) {
    // rtrim marking
    if (input_p->rtrim_nts != 0 && GET_0_8(qc_batch_p->gpu_result_p[i].counters[SIDE_NTS_QUALITY]) < input_p->min_quality)  {
        status_batch_p->read_status[i] = RTRIM_READ;
    }

    // ltrim marking
    if (input_p->ltrim_nts != 0 && GET_24_32(qc_batch_p->gpu_result_p[i].counters[SIDE_NTS_QUALITY]) < input_p->min_quality)  {
        if (status_batch_p->read_status[i] == RTRIM_READ) {
            status_batch_p->read_status[i] = TRIM_READ;
        } else {
            status_batch_p->read_status[i] = LTRIM_READ;
        }
    }
}

void filtering_read(int i, int read_length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p) {
    // i: read index within the current batch
    int length, num_n, num_out_of_quality;

    // if left or right side does not meet quality threshold the read is invalidated
    if ((GET_8_16(qc_batch_p->gpu_result_p[i].counters[SIDE_NTS_QUALITY]) < input_p->min_quality) || (GET_16_24(qc_batch_p->gpu_result_p[i].counters[SIDE_NTS_QUALITY]) < input_p->min_quality)) {
        status_batch_p->read_status[i] = INVALID_READ;
        return;
    }

    // validation variables are obtained depending on the type of trimming
    if (status_batch_p->read_status[i] == RTRIM_READ) {
        length = read_length - input_p->rtrim_nts;
        num_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[N]);
        num_out_of_quality = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[OUT_OF_QUALITY]);
    } else if (status_batch_p->read_status[i] == TRIM_READ) {
        length = read_length - input_p->rtrim_nts - input_p->ltrim_nts;
        num_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[N]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[N]);
        num_out_of_quality = GET_0_12(qc_batch_p->gpu_result_p[i].counters[OUT_OF_QUALITY]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[OUT_OF_QUALITY]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[OUT_OF_QUALITY]);
    } else if (status_batch_p->read_status[i] == LTRIM_READ) {
        length = read_length - input_p->ltrim_nts;
        num_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[N]);
        num_out_of_quality = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[OUT_OF_QUALITY]);
    } else {
        length = read_length;
        num_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]);
        num_out_of_quality = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]);
    }

    // validation filters are aplicated
    if ((length < input_p->min_read_length) ||
        (length > input_p->max_read_length) ||
        (num_n > input_p->max_n_per_read) ||
        (num_out_of_quality > input_p->max_nts_mismatch) ||
        (GET_0_8(qc_batch_p->gpu_result_p[i].counters[MEAN_MEDIAN_QUALITY]) < input_p->min_quality) ||
        (GET_8_16(qc_batch_p->gpu_result_p[i].counters[MEAN_MEDIAN_QUALITY]) < input_p->min_quality)) {

        status_batch_p->read_status[i] = INVALID_READ;
    } else {
        status_batch_p->read_status[i] = VALID_READ;
    }
}

void qc_per_read(int i, int length, int base_quality, int* mean_read_quality, long* read_sum_length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p, qc_report_t* qc_report, int cpu_num_threads) {
    int index, counters_a, counters_c, counters_g, counters_t, counters_n;
    int source_id = qc_batch_p->source_id;

    if (length > 0) {
        if (length > MAX_LINE_LENGTH) {
            char log_message[50];
            sprintf(log_message, "Be careful, length = %i > MAX_LINE_LENGTH (%i)\n", length, MAX_LINE_LENGTH);
            LOG_WARN(log_message);
        }

        if (length < qc_report[source_id].min_read_length) {
            qc_report[source_id].min_read_length = length;
        }
        if (length > qc_report[source_id].max_read_length) {
            qc_report[source_id].max_read_length = length;
        }

        read_sum_length[source_id] += length;

        // variables are filled depending on the read status
        if (status_batch_p->read_status[i] == RTRIM_READ) {
            *mean_read_quality = GET_16_24(qc_batch_p->gpu_result_p[i].counters[TRIM_QUALITY]) - base_quality;
            counters_a = GET_0_12(qc_batch_p->gpu_result_p[i].counters[A]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[A]);
            counters_c = GET_0_12(qc_batch_p->gpu_result_p[i].counters[C]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[C]);
            counters_g = GET_0_12(qc_batch_p->gpu_result_p[i].counters[G]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[G]);
            counters_t = GET_0_12(qc_batch_p->gpu_result_p[i].counters[T]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[T]);
            counters_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[N]);
        } else if (status_batch_p->read_status[i] == TRIM_READ) {
            *mean_read_quality = GET_8_16(qc_batch_p->gpu_result_p[i].counters[TRIM_QUALITY]) - base_quality;
            counters_a = GET_0_12(qc_batch_p->gpu_result_p[i].counters[A]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[A]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[A]);
            counters_c = GET_0_12(qc_batch_p->gpu_result_p[i].counters[C]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[C]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[C]);
            counters_g = GET_0_12(qc_batch_p->gpu_result_p[i].counters[G]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[G]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[G]);
            counters_t = GET_0_12(qc_batch_p->gpu_result_p[i].counters[T]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[T]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[T]);
            counters_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_12_22(qc_batch_p->gpu_result_p[i].counters[N]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[N]);
        } else if (status_batch_p->read_status[i] == LTRIM_READ) {
            *mean_read_quality = GET_24_32(qc_batch_p->gpu_result_p[i].counters[TRIM_QUALITY]) - base_quality;
            counters_a = GET_0_12(qc_batch_p->gpu_result_p[i].counters[A]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[A]);
            counters_c = GET_0_12(qc_batch_p->gpu_result_p[i].counters[C]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[C]);
            counters_g = GET_0_12(qc_batch_p->gpu_result_p[i].counters[G]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[G]);
            counters_t = GET_0_12(qc_batch_p->gpu_result_p[i].counters[T]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[T]);
            counters_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]) - GET_22_32(qc_batch_p->gpu_result_p[i].counters[N]);
        } else {
            *mean_read_quality = GET_0_8(qc_batch_p->gpu_result_p[i].counters[TRIM_QUALITY]) - base_quality;
            counters_a = GET_0_12(qc_batch_p->gpu_result_p[i].counters[A]);
            counters_c = GET_0_12(qc_batch_p->gpu_result_p[i].counters[C]);
            counters_g = GET_0_12(qc_batch_p->gpu_result_p[i].counters[G]);
            counters_t = GET_0_12(qc_batch_p->gpu_result_p[i].counters[T]);
            counters_n = GET_0_12(qc_batch_p->gpu_result_p[i].counters[N]);
        }

        qc_report[source_id].mean_read_quality += *mean_read_quality;
        qc_report[source_id].per_sequence_quality[*mean_read_quality]++;
        qc_report[source_id].a_perc += counters_a;
        qc_report[source_id].c_perc += counters_c;
        qc_report[source_id].g_perc += counters_g;
        qc_report[source_id].t_perc += counters_t;
        qc_report[source_id].n_perc += counters_n;

        qc_report[source_id].gc_histogram[100 * (counters_g + counters_c)/(counters_a + counters_c + counters_g + counters_t + counters_n)]++;

        // intervals of nt count must be calculated depending on the read status
        int length_after_trims, init_pos, end_pos;

        if (status_batch_p->read_status[i] == RTRIM_READ) {
            length_after_trims = length - input_p->rtrim_nts;
            init_pos = 0;
            end_pos = length - input_p->rtrim_nts;
        } else if (status_batch_p->read_status[i] == TRIM_READ) {
            length_after_trims = length - input_p->rtrim_nts - input_p->ltrim_nts;
            init_pos = input_p->ltrim_nts;
            end_pos = length - input_p->rtrim_nts;
        } else if (status_batch_p->read_status[i] == LTRIM_READ) {
            length_after_trims = length - input_p->ltrim_nts;
            init_pos = input_p->ltrim_nts;
            end_pos = length;
        } else {
            length_after_trims = length;
            init_pos = 0;
            end_pos = length;
        }

        //end position relative to init position
        int end_relative_pos = end_pos - init_pos;

        //offset of the quality index: init_pos + vector index + offset of the read/quality pair in data vector
        index = init_pos + qc_batch_p->read_p->data_indices[i] - 1;

        for (int j = 0 ; j < end_relative_pos ; j++) {
            index++;
            qc_report[source_id].mean_nt_quality[j] += (int) qc_batch_p->read_p->quality[index];
            qc_report[source_id].nt_counter[j]++;
        }

        // calculation of the index offset, only dependant on the order of the read [i]
        int index_offset = init_pos + qc_batch_p->read_p->data_indices[i] - 1;

        for (int j = init_pos; j < end_pos; j++) {
            index = (j * INT_COUNTERS_SIZE_IN_MEMORY) + (qc_batch_p->read_p->seq[index_offset++] & 7);

            if (status_batch_p->read_status[i] == INVALID_READ) {
                qc_batch_p->gpu_nt_type_invalid_counter_p[index]++;
            } else {
                qc_batch_p->gpu_nt_type_valid_counter_p[index - init_pos]++;
            }
        }
    }
}

void kmers_per_read(int i, int length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, status_batch_t* status_batch_p, pthread_mutex_t kmers_lock, int cpu_num_threads) {
    int k, index;
    char nt[5];
    int init_pos, end_kmers_pos;

    // intervals of nt count must be calculated depending on the read status
    if (status_batch_p->read_status[i] == RTRIM_READ) {
        init_pos = 0;
        end_kmers_pos = length - input_p->rtrim_nts - 5;
    } else if (status_batch_p->read_status[i] == TRIM_READ) {
        init_pos = input_p->ltrim_nts;
        end_kmers_pos = length - input_p->rtrim_nts - 5;
    } else if (status_batch_p->read_status[i] == LTRIM_READ) {
        init_pos = input_p->ltrim_nts;
        end_kmers_pos = length - 5;
    } else {
        init_pos = 0;
        end_kmers_pos = length - 5;
    }

    // kmers calculations in CPU
    k = qc_batch_p->read_p->data_indices[i];

    for (int d = init_pos; (d <= end_kmers_pos); d++) {
        strncpy(nt, &(qc_batch_p->read_p->seq[k+d]), 5);
        if ((nt[0] != 78) && (nt[1] != 78) && (nt[2] != 78) && (nt[3] != 78) && (nt[4] != 78)) {
            index = kmers_index(nt);

            if (status_batch_p->read_status[i] == INVALID_READ) {
                pthread_mutex_lock(&(kmers_lock));
                qc_batch_p->gpu_kmers_invalid_p[index].total_count++;
                pthread_mutex_unlock(&(kmers_lock));

                qc_batch_p->gpu_kmers_invalid_p[index].position_count[d]++;
            } else {
                pthread_mutex_lock(&(kmers_lock));
                qc_batch_p->gpu_kmers_valid_p[index].total_count++;
                pthread_mutex_unlock(&(kmers_lock));

                qc_batch_p->gpu_kmers_valid_p[index].position_count[d]++;
            }
        }
    }
}

void kmers_thread_per_read(int i, int length, results_server_input_t* input_p, qc_batch_t* qc_batch_p, int* kmers_valid_thread_counter, int* kmers_invalid_thread_counter, status_batch_t* status_batch_p, int cpu_num_threads) {
    int k, index;
    char nt[5];
    int init_pos, end_kmers_pos;

    // intervals of nt count must be calculated depending on the read status
    if (status_batch_p->read_status[i] == RTRIM_READ) {
        init_pos = 0;
        end_kmers_pos = length - input_p->rtrim_nts - 5;
    } else if (status_batch_p->read_status[i] == TRIM_READ) {
        init_pos = input_p->ltrim_nts;
        end_kmers_pos = length - input_p->rtrim_nts - 5;
    } else if (status_batch_p->read_status[i] == LTRIM_READ) {
        init_pos = input_p->ltrim_nts;
        end_kmers_pos = length - 5;
    } else {
        init_pos = 0;
        end_kmers_pos = length - 5;
    }

    // kmers calculations in CPU
    k = qc_batch_p->read_p->data_indices[i];

    for (int d = init_pos; (d <= end_kmers_pos); d++) {
        strncpy(nt, &(qc_batch_p->read_p->seq[k+d]), 5);
        if ((nt[0] != 78) && (nt[1] != 78) && (nt[2] != 78) && (nt[3] != 78) && (nt[4] != 78)) {
            index = kmers_index(nt);

            if (status_batch_p->read_status[i] == INVALID_READ) {
                kmers_invalid_thread_counter[index]++;
                qc_batch_p->gpu_kmers_invalid_p[index].position_count[d]++;
            } else {
                kmers_valid_thread_counter[index]++;
                qc_batch_p->gpu_kmers_valid_p[index].position_count[d]++;
            }
        }
    }
}

void qc_per_batch(int source_id, results_server_input_t* input_p, qc_batch_t* qc_batch_p, qc_report_t* qc_report, int valid) {
    if (valid) { //valid reads
        for (int j = 0 ; j < qc_report[source_id].max_read_length ; j++) {
            qc_report[source_id].nt_type_counter[j][A] += qc_batch_p->gpu_nt_type_valid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + A];
            qc_report[source_id].nt_type_counter[j][C] += qc_batch_p->gpu_nt_type_valid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + C];
            qc_report[source_id].nt_type_counter[j][G] += qc_batch_p->gpu_nt_type_valid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + G];
            qc_report[source_id].nt_type_counter[j][T] += qc_batch_p->gpu_nt_type_valid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + T];
            qc_report[source_id].nt_type_counter[j][N] += qc_batch_p->gpu_nt_type_valid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + N];
        }

        if (input_p->kmers_on) {
            for (int j = 0; j < KMERS_COMBINATIONS; j++) {

                qc_report[source_id].qc_kmers[j].total_count += qc_batch_p->gpu_kmers_valid_p[j].total_count;

                for (int d = 0; d < MAX_LINE_LENGTH; d++) {
                    qc_report[source_id].qc_kmers[j].position_count[d] += qc_batch_p->gpu_kmers_valid_p[j].position_count[d];
                }
            }
        }
    } else {  //invalid reads
        for (int j = 0 ; j < qc_report[source_id].max_read_length ; j++) {
            qc_report[source_id].nt_type_counter[j][A] += qc_batch_p->gpu_nt_type_invalid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + A];
            qc_report[source_id].nt_type_counter[j][C] += qc_batch_p->gpu_nt_type_invalid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + C];
            qc_report[source_id].nt_type_counter[j][G] += qc_batch_p->gpu_nt_type_invalid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + G];
            qc_report[source_id].nt_type_counter[j][T] += qc_batch_p->gpu_nt_type_invalid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + T];
            qc_report[source_id].nt_type_counter[j][N] += qc_batch_p->gpu_nt_type_invalid_counter_p[(j * COUNTERS_SIZE * sizeof(int)) + N];
        }

        if (input_p->kmers_on) {
            for (int j = 0; j < KMERS_COMBINATIONS; j++) {

                qc_report[source_id].qc_kmers[j].total_count += qc_batch_p->gpu_kmers_invalid_p[j].total_count;

                for (int d = 0; d < MAX_LINE_LENGTH; d++) {
                    qc_report[source_id].qc_kmers[j].position_count[d] += qc_batch_p->gpu_kmers_invalid_p[j].position_count[d];
                }
            }
        }
    }
}

void qc_results(results_server_input_t* input_p, long* read_sum_length, qc_report_t* qc_report, int valid) {
    int i = 0;

    for (int source_id = 0; source_id < input_p->num_sources; source_id++) {
        if (time_flag) {
            start_timer(t1_result);
        }

        // calculate mean read length
        if (qc_report[source_id].nb_reads != 0) {
            qc_report[source_id].mean_read_length = read_sum_length[source_id] / qc_report[source_id].nb_reads;
        } else {
            qc_report[source_id].mean_read_length = 0;
        }

        // calculate mean quality per read
        if (qc_report[source_id].nb_reads != 0) {
            qc_report[source_id].mean_read_quality /= qc_report[source_id].nb_reads;
        } else {
            qc_report[source_id].mean_read_quality = 0;
        }

        // calculate mean quality per nucleotide position
        for (int j = 0; j < MAX_LINE_LENGTH ; j++) {
            if (qc_report[source_id].nt_counter[j] > 0) {
                qc_report[source_id].mean_nt_quality[j] /= qc_report[source_id].nt_counter[j];
                qc_report[source_id].mean_nt_quality[j] -= input_p->base_quality;
            }
        }

        // recalculate mean quality per read
        i = 0;
        for (int j = 0; j < MAX_LINE_LENGTH ; j++) {
            if (qc_report[source_id].mean_nt_quality[j] >= 0) {
                qc_report[source_id].mean_read_quality_from_nt += qc_report[source_id].mean_nt_quality[j];
                i++;
            }
        }
        qc_report[source_id].mean_read_quality_from_nt /= i;

        // calculate percentage of A, C, G, T, N nucleotides
        long nb_nt[MAX_NUM_PRODUCERS];
        nb_nt[source_id] = qc_report[source_id].a_perc + qc_report[source_id].c_perc + qc_report[source_id].g_perc + qc_report[source_id].t_perc + qc_report[source_id].n_perc;

        qc_report[source_id].a_perc = 100 * (qc_report[source_id].a_perc / nb_nt[source_id]);
        qc_report[source_id].c_perc = 100 * (qc_report[source_id].c_perc / nb_nt[source_id]);
        qc_report[source_id].g_perc = 100 * (qc_report[source_id].g_perc / nb_nt[source_id]);
        qc_report[source_id].t_perc = 100 * (qc_report[source_id].t_perc / nb_nt[source_id]);
        qc_report[source_id].n_perc = 100 * (qc_report[source_id].n_perc / nb_nt[source_id]);

        mean_reads_per_batch += qc_report[source_id].nb_reads;

        // sorting kmers by total count
        if (input_p->kmers_on) {
            qsort(qc_report[source_id].qc_kmers, KMERS_COMBINATIONS, sizeof(qc_kmers_t), kmers_sort);
        }

        if (time_flag) {
            stop_timer(t1_result, t2_result, result_time);
        }

        // and finally, print qc report, data files and graphs
        if (time_flag) {
            start_timer(t1_reporting);
        }
        generate_report(qc_report[source_id], input_p->source_p[source_id].filename, input_p->base_quality, input_p->kmers_on, input_p->cg_flag, input_p->report_directory, valid);
        if (time_flag) {
            stop_timer(t1_reporting, t2_reporting, reporting_time);
        }
    }
}

int kmers_index(char* nt) {
    int index = 0;
    index += ((nt[4] >> 1) & 3);
    index += (((nt[3] >> 1) & 3) << 2);
    index += (((nt[2] >> 1) & 3) << 4);
    index += (((nt[1] >> 1) & 3) << 6);
    index += (((nt[0] >> 1) & 3) << 8);

    return index;
}
