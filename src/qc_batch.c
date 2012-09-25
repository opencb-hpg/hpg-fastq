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


#include <stdlib.h>
#include <pthread.h>
#include <string.h>

#include "qc_batch.h"

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

void qc_batch_free(qc_batch_t* qc_batch_p, int all) {
    if (qc_batch_p == NULL) return;

    if (all) {
        if (qc_batch_p->read_p != NULL) fastq_batch_free(qc_batch_p->read_p);
        qc_batch_p->read_p = NULL;
    }
    
    if (qc_batch_p->gpu_result_p != NULL) free(qc_batch_p->gpu_result_p);
    qc_batch_p->gpu_result_p = NULL;
    if (qc_batch_p->gpu_kmers_valid_p != NULL) free(qc_batch_p->gpu_kmers_valid_p);
    qc_batch_p->gpu_kmers_valid_p = NULL;
    if (qc_batch_p->gpu_kmers_invalid_p != NULL) free(qc_batch_p->gpu_kmers_invalid_p);
    qc_batch_p->gpu_kmers_invalid_p = NULL;    
    if (qc_batch_p->gpu_nt_type_valid_counter_p != NULL) free(qc_batch_p->gpu_nt_type_valid_counter_p);
    qc_batch_p->gpu_nt_type_valid_counter_p = NULL;
    if (qc_batch_p->gpu_nt_type_invalid_counter_p != NULL) free(qc_batch_p->gpu_nt_type_invalid_counter_p);
    qc_batch_p->gpu_nt_type_invalid_counter_p = NULL;
    
    free(qc_batch_p);
    qc_batch_p = NULL;
}

void qc_kmers_init(qc_kmers_t* qc_kmers_p) {
    for (int i = 0; i < KMERS_COMBINATIONS; i++) {
        qc_kmers_p[i].id = i;
        qc_kmers_p[i].total_count = 0;
        memset(qc_kmers_p[i].position_count, 0, MAX_LINE_LENGTH * sizeof(int));
    }
}

char* kmers_string(int index, char* kmer) {
    int rest, quotient;

    quotient = index;
    kmer[5] = '\0';

    for (int i = 0; i < 5; i++) {
        rest = quotient % 4;
        quotient /= 4;

        switch (rest) {
            case 0:
                kmer[4-i] = 'A';
                break;
            case 1:
                kmer[4-i] = 'C';
                break;
            case 2:
                kmer[4-i] = 'T';
                break;
            case 3:
                kmer[4-i] = 'G';
                break;
        }
    }

    return kmer;
}

int kmers_sort(const void* k1, const void* k2) {
    qc_kmers_t* kmer1_p = (qc_kmers_t*) k1;
    qc_kmers_t* kmer2_p = (qc_kmers_t*) k2;

    return (kmer2_p->total_count - kmer1_p->total_count);
}
