/*
 * stats_fastq.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "stats_fastq.h"

//--------------------------------------------------------------------
// stats counters
//--------------------------------------------------------------------

stats_counters_t *stats_counters_new() {
  stats_counters_t *sc = calloc(1, sizeof(stats_counters_t));

  sc->num_reads = 0;
  sc->phred = QUALITY_PHRED33_VALUE;

  sc->min_length = 100000;
  sc->max_length = 0;
  sc->acc_length = 0;
  sc->mean_length = 0.0f;

  sc->acc_quality = 0;
  sc->mean_quality = 0.0f;

  sc->kh_length_histogram = kh_init(32);
  sc->kh_quality_histogram = kh_init(32);
  sc->kh_gc_histogram = kh_init(32);

  sc->kh_count_quality_per_nt = kh_init(32);
  sc->kh_acc_quality_per_nt = kh_init(32);
  
  sc->kh_num_As_per_nt = kh_init(32);
  sc->kh_num_Cs_per_nt = kh_init(32);
  sc->kh_num_Ts_per_nt = kh_init(32);
  sc->kh_num_Gs_per_nt = kh_init(32);
  sc->kh_num_Ns_per_nt = kh_init(32);

  return sc;
}

void stats_counters_free(stats_counters_t *sc) {
  if (sc) {
    
    kh_destroy(32, sc->kh_length_histogram);
    kh_destroy(32, sc->kh_quality_histogram);
    kh_destroy(32, sc->kh_gc_histogram);

    kh_destroy(32, sc->kh_count_quality_per_nt);
    kh_destroy(32, sc->kh_acc_quality_per_nt);
    
    kh_destroy(32, sc->kh_num_As_per_nt);
    kh_destroy(32, sc->kh_num_Cs_per_nt);
    kh_destroy(32, sc->kh_num_Ts_per_nt);
    kh_destroy(32, sc->kh_num_Gs_per_nt);
    kh_destroy(32, sc->kh_num_Ns_per_nt);
    
    free(sc);
  }
}

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct fastq_stats_wf_batch {
  stats_options_t *options;
  fastq_read_stats_options_t *read_stats_options;
  array_list_t *fq_reads;
  array_list_t *fq_stats;
  stats_counters_t *stats_counters;
} fastq_stats_wf_batch_t;

fastq_stats_wf_batch_t *fastq_stats_wf_batch_new(stats_options_t *opts,
						 fastq_read_stats_options_t *read_stats_options,
						 array_list_t *fq_reads,
						 array_list_t *fq_stats,
						 stats_counters_t *stats_counters) {
  
  fastq_stats_wf_batch_t *b = (fastq_stats_wf_batch_t *) calloc(1, 
								sizeof(fastq_stats_wf_batch_t));
  
  b->options = opts;
  b->read_stats_options = read_stats_options;
  b->fq_reads = fq_reads;
  b->fq_stats = fq_stats;
  b->stats_counters = stats_counters;
  
  return b;
}

void fastq_stats_wf_batch_free(fastq_stats_wf_batch_t *b) {
  if (b) free(b);
}

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct fastq_stats_wf_input {
  stats_options_t *options;
  fastq_file_t *in_file;
  fastq_read_stats_options_t *read_stats_options;
  stats_counters_t *stats_counters;
} fastq_stats_wf_input_t;

fastq_stats_wf_input_t *fastq_stats_wf_input_new(stats_options_t *opts,
						 fastq_file_t *in_file,
						 stats_counters_t *stats_counters) {
  
  fastq_stats_wf_input_t *wfi = (fastq_stats_wf_input_t *) calloc(1, sizeof(fastq_stats_wf_input_t));

  wfi->in_file = in_file;
  wfi->options = opts;
  wfi->stats_counters = stats_counters;
  wfi->read_stats_options = fastq_read_stats_options_new(0, 10000, 1);

  return wfi;
}

void fastq_stats_wf_input_free(fastq_stats_wf_input_t *wfi) {
  if (wfi) {

    if (wfi->read_stats_options) fastq_read_stats_options_free(wfi->read_stats_options);

    free(wfi);
  }
}

//--------------------------------------------------------------------
// workflow producer : read FastQ file
//--------------------------------------------------------------------

void *fastq_stats_producer(void *input) {
  
  fastq_stats_wf_input_t *wf_input = (fastq_stats_wf_input_t *) input;
  fastq_stats_wf_batch_t *wf_batch = NULL;
  size_t max_num_reads = wf_input->options->batch_size;

  
  array_list_t *fq_reads = array_list_new(max_num_reads,
					  1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  size_t num_reads = fastq_fread_se(fq_reads, max_num_reads, wf_input->in_file);

  if (num_reads) {
    array_list_t *fq_stats = array_list_new(num_reads,
					    1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

    wf_batch = fastq_stats_wf_batch_new(wf_input->options,
					wf_input->read_stats_options,
					fq_reads,
					fq_stats,
					wf_input->stats_counters);				      
  } else {
    array_list_free(fq_reads, NULL);
  }
    
  return wf_batch;
}

//--------------------------------------------------------------------
// workflow worker : compute statistics per read
//--------------------------------------------------------------------

int fastq_stats_worker(void *data) {
  fastq_stats_wf_batch_t *batch = (fastq_stats_wf_batch_t *) data;

  fastq_reads_stats(batch->fq_reads,
		    batch->read_stats_options,
		    batch->fq_stats);

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : merge statistics
//--------------------------------------------------------------------

int fastq_stats_consumer(void *data) {
  fastq_stats_wf_batch_t *batch = (fastq_stats_wf_batch_t *) data;

  khiter_t k;
  khash_t(32) *h;

  fastq_read_t *fq_read;
  fastq_read_stats_t *read_stats;
  stats_counters_t *counters = batch->stats_counters;

  int ret;
  size_t value, read_length;
  size_t num_reads = array_list_size(batch->fq_stats);
  
  for (size_t i = 0 ; i < num_reads; i++) {
    fq_read = array_list_get(i, batch->fq_reads);
    read_stats = array_list_get(i, batch->fq_stats);
    if (read_stats) {
      counters->num_reads++;
      
      // length, min, max, acc
      read_length = read_stats->length;
      counters->acc_length += read_length;
      if (counters->min_length > read_length) 
	counters->min_length = read_length;
      if (counters->max_length < read_length) 
	counters->max_length = read_length;

      // quality
      counters->acc_quality += read_stats->quality_average;
      
      // A, C, T, G, N
      counters->num_As += read_stats->num_A;
      counters->num_Cs += read_stats->num_C;
      counters->num_Ts += read_stats->num_T;
      counters->num_Gs += read_stats->num_G;
      counters->num_Ns += read_stats->Ns;

      // length histogram
      value = read_length;
      h = counters->kh_length_histogram;
      k = kh_put(32, h, value, &ret);
      if (ret) {
	kh_value(h, k) = 1;
      } else {
	kh_value(h, k) =  kh_value(h, k) + 1;
      }

      // quality histogram
      value = round(read_stats->quality_average);
      h = counters->kh_quality_histogram;
      k = kh_put(32, h, value, &ret);
      if (ret) {
	kh_value(h, k) = 1;
      } else {
	kh_value(h, k) =  kh_value(h, k) + 1;
      }

      // GC histogram
      value = 100 * (read_stats->num_G + read_stats->num_C) / read_length;
      h = counters->kh_gc_histogram;
      k = kh_put(32, h, value, &ret);
      if (ret) {
	kh_value(h, k) = 1;
      } else {
	kh_value(h, k) =  kh_value(h, k) + 1;
      }


      // stats per nt position
      for (size_t j = 0; j < read_length; j++) {

	// quality counter
	h = counters->kh_count_quality_per_nt;
	k = kh_put(32, h, j, &ret);
	if (ret) {
	  kh_value(h, k) = 1;
	} else {
	  kh_value(h, k) =  kh_value(h, k) + 1;
	}

	// quality accumulator 
	h = counters->kh_acc_quality_per_nt;
	k = kh_put(32, h, j, &ret);
	if (ret) {
	  kh_value(h, k) = fq_read->quality[j];
	} else {
	  kh_value(h, k) =  kh_value(h, k) + fq_read->quality[j];
	}

	// A, T, C, G, N
	h = NULL;
	switch (fq_read->sequence[j]) {
	case 'A': 	h = counters->kh_num_As_per_nt;
	  break;
	case 'T':	h = counters->kh_num_Ts_per_nt;
	  break;
	case 'C':	h = counters->kh_num_Cs_per_nt;
	  break;
	case 'G':	h = counters->kh_num_Gs_per_nt;
	  break;
	case 'N':	h = counters->kh_num_Ns_per_nt;
	  break;
	default:	break;
	}

	if (h) {
	  k = kh_put(32, h, j, &ret);
	  if (ret) {
	    kh_value(h, k) = 1;
	  } else {
	    kh_value(h, k) =  kh_value(h, k) + 1;
	  }
	}
      }
    }
  }
  
  // free memory
  fastq_stats_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void stats_fastq(stats_options_t *opts) {

  fastq_file_t *fq_file = fastq_fopen(opts->in_filename);
  stats_counters_t *counters = stats_counters_new();
  counters->phred = opts->quality_encoding_value;
  
  //------------------------------------------------------------------
  // workflow management
  //
  fastq_stats_wf_input_t *wf_input = fastq_stats_wf_input_new(opts,
							      fq_file,
							      counters);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {fastq_stats_worker};
  char *stage_labels[] = {"FASTQ stats worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(fastq_stats_producer, "FASTQ stats producer", wf);
  workflow_set_consumer(fastq_stats_consumer, "FASTQ stats consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);
  
  // free memory
  workflow_free(wf);
  fastq_stats_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  stats_report(opts->in_filename, counters, opts->out_dirname);

  // free memory
  stats_counters_free(counters);
  fastq_fclose(fq_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
