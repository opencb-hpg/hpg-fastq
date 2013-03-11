/*
 * main.c
 *
 *  Created on: Mar 8, 2013
 *      Author: jtarraga
 */

//#include "bioformats/bam-sam/fast_stats.h"
//#include "bioformats/bam-sam/bam_stats_report.h"
//#include "bioformats/bam-sam/bam_filter.h"

//#include "sqlite3.h"

#include "stats_options.h"
//#include "filter_options.h"
//#include "index_options.h"
//#include "sort_options.h"

//------------------------------------------------------------------------

void usage(char *exec_name);

extern int bam_index_build2(const char *fn, const char *_fnidx);

//------------------------------------------------------------------------
//                    M A I N     F U N C T I O N
//------------------------------------------------------------------------

int main (int argc, char *argv[]) {
  char *exec_name = argv[0];  
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    usage(exec_name);
  }

  char *command_name = argv[1];  

  argc--;
  argv++;

  if (strcmp(command_name, "stats") == 0) {

    //--------------------------------------------------------------------
    //                  S T A T S     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display stats options
    stats_options_t *opts = parse_stats_options(exec_name, command_name, 
						argc, argv);
    validate_stats_options(opts);
    display_stats_options(opts);
    
    // run and display stats
    stats_fastq(opts);

    // free memory
    free_stats_options(opts);

  } else if (strcmp(command_name, "filter" ) == 0) {

    //--------------------------------------------------------------------
    //                  F I L T E R     C O M M A N D
    //--------------------------------------------------------------------

    printf("Not implemented yet !!\n");
    /*
    // parse, validate and display filter options
    filter_options_t *opts = parse_filter_options(exec_name, command_name, 
						  argc, argv);
    validate_filter_options(opts);
    display_filter_options(opts);
    
    // run filter
    filter_bam(opts);

    // free memory
    free_filter_options(opts);
    */
  } else {

    //--------------------------------------------------------------------
    //                  U N K N O W N     C O M M A N D
    //--------------------------------------------------------------------

    usage(exec_name);
  }
  printf("Done !\n");
}

//------------------------------------------------------------------------

void usage(char *exec_name) {
    printf("Program: %s (High-performance tools for handling BAM files)\n", exec_name);
    printf("Version: 1.0.0\n");
    printf("\n");
    printf("Usage: %s <command> [options]\n", exec_name);
    printf("\n");
    printf("Command: stats\t\tstatistics summary\n");
    printf("         filter\t\tfilter a FastQ file by using advanced criteria\n");
    printf("\n");    
    printf("For more information about a certain command, type %s <command> --help\n", exec_name);
    exit(-1);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
