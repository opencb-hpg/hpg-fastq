#ifndef STATS_OPTIONS_H
#define STATS_OPTIONS_H

/*
 * stats_options.h
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

//#define NUM_STATS_OPTIONS	9
#define NUM_STATS_OPTIONS	7

//------------------------------------------------------------------------

typedef struct stats_options { 
  int log_level;
  int verbose;
  int help;
  int num_threads;
  int batch_size;

  char* in_filename;
  char* out_dirname;
  char* gff_region_filename;
  char* region_list;

  char *exec_name;
  char *command_name;
} stats_options_t;

//------------------------------------------------------------------------

stats_options_t *new_stats_options(char *exec_name, char *command_nane);

stats_options_t *parse_stats_options(char *exec_name, char *command_nane,
				     int argc, char **argv);

void free_stats_options(stats_options_t *opts);

void validate_stats_options(stats_options_t *opts);

void display_stats_options(stats_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
