
#include <stdlib.h>
#include <string.h>

#include "prepro_batch.h"

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

void status_batch_list_init(status_batch_list_t* list_p) {
    memset(list_p, 0, sizeof(status_batch_list_t));
    list_p->length = 0;
    list_p->first_p = NULL;
    list_p->last_p = NULL;
    pthread_mutex_init(&(list_p->lock), NULL);
}

void status_batch_list_insert(status_batch_t *status_batch_p, status_batch_list_t *status_batch_list_p) {
    if (status_batch_list_p == NULL) return;

    pthread_mutex_lock(&status_batch_list_p->lock);

    if (status_batch_list_p->first_p == NULL) {
        status_batch_p->prev_p = NULL;
        status_batch_p->next_p = NULL;
        status_batch_list_p->first_p = status_batch_p;
        status_batch_list_p->last_p = status_batch_p;
    } else {
        status_batch_list_p->last_p->next_p = status_batch_p;
        status_batch_p->prev_p = status_batch_list_p->last_p;
        status_batch_p->next_p = NULL;
        status_batch_list_p->last_p = status_batch_p;
    }
    status_batch_list_p->length++;
    status_batch_list_p->length_by_source_id[status_batch_p->source_id]++;

    pthread_mutex_unlock(&status_batch_list_p->lock);
}

status_batch_t* status_batch_list_remove(status_batch_list_t* status_batch_list_p) {
    if (status_batch_list_p == NULL) return NULL;

    pthread_mutex_lock(&status_batch_list_p->lock);

    // just get the first element, and if is not null update the first-element pointer
    status_batch_t* status_batch_p = status_batch_list_p->first_p;
    if (status_batch_p != NULL) {
        status_batch_list_p->first_p = status_batch_p->next_p;
        status_batch_list_p->length--;
    }

    pthread_mutex_unlock(&status_batch_list_p->lock);

    return status_batch_p;
}

status_batch_t* status_batch_list_get_next_by_id(status_batch_list_t* list_p, int source_id, int batch_id) {
    if (list_p == NULL) {
        return NULL;
    }

    pthread_mutex_lock(&list_p->lock);

    status_batch_t* batch_p = list_p->first_p;

    while (batch_p != NULL) {
        if (batch_p->source_id == source_id && batch_p->id == batch_id) {
            if (list_p->first_p == batch_p && list_p->last_p == batch_p) {
                list_p->first_p = NULL;
                list_p->last_p = NULL;
            } else if (list_p->first_p == batch_p) {
                list_p->first_p = batch_p->next_p;
                list_p->first_p->prev_p = NULL;
            } else if (list_p->last_p == batch_p) {
                list_p->last_p = batch_p->prev_p;
                list_p->last_p->next_p = NULL;
            } else {
                batch_p->prev_p->next_p = batch_p->next_p;
                batch_p->next_p->prev_p = batch_p->prev_p;
            }

            list_p->length--;
            (list_p->length_by_source_id[source_id])--;
            break;
        }

        batch_p = batch_p->next_p;
    }

    pthread_mutex_unlock(&list_p->lock);

    return batch_p;
}

void status_batch_list_print(status_batch_list_t* list_p) {
    int nt_count;

    pthread_mutex_lock(&list_p->lock);

    printf("Number of items: %i\n", list_p->length);

    status_batch_t* status_batch_p = list_p->first_p;

    int i = 0;
    while (status_batch_p != NULL) {
        printf("+++ %i, status_batch_source_id: %i, status_batch id: %i, num_reads: %i\n", i++, status_batch_p->source_id, status_batch_p->id, status_batch_p->num_reads);
        status_batch_p = status_batch_p->next_p;
    }

    pthread_mutex_unlock(&list_p->lock);
}

void status_batch_list_items_free(status_batch_list_t* list_p) {
    status_batch_t* item_p;

    pthread_mutex_lock(&list_p->lock);

    while ((item_p = status_batch_list_remove(list_p)) != NULL) {
        list_p->length--;
        status_batch_free(item_p, 1);
    }

    pthread_mutex_unlock(&list_p->lock);
}

void status_batch_free(status_batch_t* status_batch_p, int all) {
    if (status_batch_p == NULL) return;

    if (all) {
        if (status_batch_p->read_status != NULL) free(status_batch_p->read_status);
    }
    free(status_batch_p);
}

int calculate_last_nt_mean_quality(char* quality, int* data_indices, int position, int last_nts) {
    int i, acc_quality = 0;

    for (i = 0; i <= last_nts; i++) {
        acc_quality += quality[data_indices[position + 1] - i - 1];
    }

    return acc_quality / last_nts;
}
