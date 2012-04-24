
#ifndef QC_REPORT_H
#define QC_REPORT_H

#include <stdio.h>

#include "chaos_game.h"
#include "qc_batch.h"

#define QC_SUFFIX		".qc"
#define CG_SUFFIX		".cg"
#define POSITION_DATA_SUFFIX    ".pos.dat"
#define QUALITY_DATA_SUFFIX     ".quality.dat"
#define GC_DATA_SUFFIX          ".gc.dat"
#define KMERS_TOTAL_SUFFIX	".kmers.dat"
#define KMERS_POSITION_SUFFIX	".kmers.pos.dat"
#define HTML_FILE_SUFFIX	".qc.html"

#define POS_COLUMN	        1
#define NT_QUALITY_COLUMN	2
#define NUM_READ_LENGTH_COLUMN	3 //This column contains the number of reads with a given length for read length histogram
#define A_COLUMN		4
#define C_COLUMN		5
#define G_COLUMN		6
#define T_COLUMN		7
#define N_COLUMN		8
#define GC_COLUMN		9

#define READ_QUALITY_COLUMN	2
#define GC_HISTOGRAM_COLUMN	2

#define KMER_POS_COLUMN		1
#define KMER_1_COLUMN		2
#define KMER_2_COLUMN		3
#define KMER_3_COLUMN		4
#define KMER_4_COLUMN		5
#define KMER_5_COLUMN		6

#define GRAPH_MAX_COLUMNS	10
#define MAX_TITLE_LENGTH	25
#define KMERS_FILE_BUFFER_SIZE	227

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Structure for qc report data 
* 
* Structure containing data for the final QC report. One structure for file and one for valid and invalid reads.
*/
typedef struct qc_report {
    unsigned long nb_reads;				/**< Number of reads in the fastq reads processed. */
    unsigned int min_read_length;			/**< Minimum length of the processed reads. */
    unsigned int max_read_length;			/**< Maximum length of the processed reads. */
    unsigned int mean_read_length;			/**< Mean read length. */
    unsigned long mean_read_quality;			/**< Mean read quality (calculated read by read). */
    int mean_read_quality_from_nt;			/**< Mean read quality (calculated nt by nt). */
    int nt_counter[MAX_LINE_LENGTH];			/**< Nt counter by nucleotide position. */
    long mean_nt_quality[MAX_LINE_LENGTH];		/**< Mean quality by nucleotide position. */
    int per_sequence_quality[MAX_QUALITY_VALUE];	/**< Histogram of read qualities. */
    int nt_type_counter[MAX_LINE_LENGTH][COUNTERS_SIZE];/**< Minimum length of the precessed reads. */
    int gc_histogram[101];				/**< Histogram of % of GC bases in the read. */
    qc_kmers_t* qc_kmers;				/**< Pointer to qc kmers structure. */
    float a_perc;					/**< %A. */ 
    float c_perc;					/**< %C. */ 
    float g_perc;					/**< %G. */ 
    float t_perc;					/**< %T. */
    float n_perc;					/**< %N. */
  
    //Fields for chaos game
    char pgm_fastq_filename[MAX_PGM_FILENAME_LENGTH];	/**< Filename of the fastq sequence file signature. */	
    char pgm_quality_filename[MAX_PGM_FILENAME_LENGTH];	/**< Filename of the fastq quality file signature. */
    char pgm_diff_filename[MAX_PGM_FILENAME_LENGTH];	/**< Filename of the normalized diff between fastq file and reference genome signature. */
    int dim_n;						/**< Dimension of the tables (dim_n X dim_n). */
    int word_size_k;					/**< Word size. */
    unsigned int fq_word_count;				/**< Total read words. */
    double mean_table_dif_value;			/**< Mean value of the diff matrix. */
    double standard_deviation_table_dif_value;		/**< Standard deviation value of the diff matrix. */
    int highest_table_dif_value;			/**< Highest value in the diff matrix. */
    int lowest_table_dif_value;				/**< Lowest value in the diff matrix. */
} qc_report_t;

/**
* @brief Structure for QC graphs parameters
* 
* Structure containing parameters for graphics representation of QC graphs
*/
typedef struct qc_graph {
    int x_autoscale;					/**< Autoscale flag for X axis. */
    int x_start;					/**< X axis start coordinate. */
    int x_end;						/**< X axis end coordinate. */
    int y_autoscale;					/**< Autoscale flag for Y axis. */
    int y_start;					/**< Y axis start coordinate. */
    int y_end;						/**< Y axis end coordinate. */
    int lmargin;					/**< Left margin. */
    int rmargin;					/**< Right margin. */	
    int tmargin;					/**< Top margin. */
    int bmargin;					/**< Bottom margin. */
    int num_y_columns;					/**< Number of columns in the Y axis. */
    int x_column;					/**< Number of column in data file with X axis data. */
    int y_columns[GRAPH_MAX_COLUMNS];			/**< Numbers of columns in data file that are represented in Y axis. */
    char* y_titles[GRAPH_MAX_COLUMNS];			/**< Titles of series represented in Y axis. */
    char* title;					/**< Title of the graphic. */
    char* xlabel;					/**< Label of X axis. */	
    char* ylabel;					/**< Label of Y axis. */
    char* type; 					/**< Type of graph. */
} qc_graph_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Fills qc_report from chaos_game_data structure
*  @param chaos_game_data_p pointer to the chaos game data structure
*  @param[in,out] qc_report_p pointer to the qc_report structure
*  @return void
*  
*  Fills proper qc_report fields from chaos_game_data structure
*/
void qc_report_get_values_from_chaos_game(chaos_game_data_t* chaos_game_data_p, qc_report_t* qc_report);

/**
*  @brief Generates HTML and text reports and data files
*  @param qc_report pointer to the qc_report structure
*  @param inputfilename filename of the processed fastq file 
*  @param base_quality base quality for quality normalization
*  @param kmers_on flag for kmers (0: kmers off, 1: kmers on)
*  @param cg_flag chaos game flag (0: no chaos game, 1: chaos game executed)
*  @param report_directory directory where output files will be written
*  @param valid flag of report for valid or invalid reads (0: invalid, 1: valid)
*  @return void
*  
*  Generates the output of QC process: HTML and text reports and data files
*/
void generate_report(qc_report_t qc_report, char* inputfilename, int base_quality, int kmers_on, int cg_flag, char* report_directory, int valid);

//void plot_nt_quality_public(char* data_filename, char* graph_filename);

#endif	/*  QC_REPORT_H	  */















