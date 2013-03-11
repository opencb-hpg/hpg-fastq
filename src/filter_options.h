#ifndef FILTER_OPTIONS_H
#define FILTER_OPTIONS_H

/*
 * filter_options.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "argtable2.h"
#include "libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//============================ DEFAULT VALUES ============================

//------------------------------------------------------------------------

//#define NUM_FILTER_OPTIONS     16
#define NUM_FILTER_OPTIONS     14

//------------------------------------------------------------------------

typedef struct filter_options {
  int log_level;
  int verbose;
  int help;
  int num_threads;
  int batch_size;

  int mapped;
  int unmapped;
  int unique;
  int proper_pairs;
  int min_num_errors;
  int max_num_errors;
  int min_quality;
  int max_quality;
  int min_length;
  int max_length;

  char *length_range;
  char *quality_range;
  char *num_errors_range;

  char* in_filename;
  char* out_dirname;
  char* gff_region_filename;
  char* region_list;

  char *exec_name;
  char *command_name;
} filter_options_t;

//------------------------------------------------------------------------

filter_options_t *new_filter_options(char *exec_name, char *command_nane);

filter_options_t *parse_filter_options(char *exec_name, char *command_nane,
				       int argc, char **argv);

void free_filter_options(filter_options_t *opts);

void validate_filter_options(filter_options_t *opts);

void display_filter_options(filter_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
