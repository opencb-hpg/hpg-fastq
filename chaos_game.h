
#ifndef CHAOS_GAME_H
#define CHAOS_GAME_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "qc_batch.h"
#include "prepro.h"
#include "prepro_batch.h"
#include "fastq_commons.h"

#define NUM_LETTERS           		4  		//Number of letters used for base encoding
#define N_UNDEFINED_BASE  		5		//Associated constant to N 
#define GREATER_THAN           		62		//Associated constant to > symbol
#define MAX_QUALITY_IN_TABLE   		62  		//Maximum quality value of the PHRED SCORES
#define EPSILON     			0.00001		//Computers do not have infinite real numbers so it is necessary 
							//to define an infinitesimal discrete EPSILON value as follows:
							// 1/(2 exp 14) = 0.00006

#define MIN_K_IMAGE_VALUE		7		//if k < MIN image is resized to make it more accesible to human sight
#define MIN_IMAGE_PIXEL_SIZE		128		//Image size: 128 X 128

#define MAX_PGM_FILENAME_LENGTH		100		//Size to reserve for .pgm filenames
#define	K_VALUE_INFIX			"_k="		//Infix to append for k value 
#define	FASTQ_PGM_FILENAME_SUFFIX	"_FG.pgm"	//Suffix and extension to add in fastq pgm file
#define	QUALITY_PGM_FILENAME_SUFFIX	"_QQ.pgm"	//Suffix and extension to add in qualities pgm file
#define	DIFF_PGM_FILENAME_SUFFIX	"_FG_dif.pgm"	//Suffix and extension to add in diff pgm file

#define uchar  unsigned char

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Header of Genomic Signature
* 
* Structure containing the fields of the header of a Genomic Signature file
*/
typedef struct header_gs {
    char          gs_filename[180];	/**< Genomic Signature filename. */
    unsigned int  word_size_k;		/**< Word size (k). */
    unsigned int  dim_x, dim_y;		/**< X and Y dimensions. */
    unsigned int  ref_word_count; 	/**< Total words of the reference sequence. */
} header_gs_t;

/**
* @brief Chaos Game data
* 
* Structure containing data used in the Chaos Game (filenames, matrix, work size, ...)
*/
typedef struct chaos_game_data {
    unsigned int**  table_seq;				/**< Table for sequences. */
    unsigned int**  table_q;				/**< Table for qualities. */
    unsigned int**  table_gs;				/**< Table for genomic signature. */
    int**           table_dif; 				/**< Table for diff between TABLE_SEQ and TABLE_GS. */
    double 	  mean_table_dif_value;			/**< Mean value of the diff matrix. */
    double 	  standard_deviation_table_dif_value;	/**< Standard deviation value of the diff matrix. */
    int 		  highest_table_dif_value;	/**< Highest value in the diff matrix. */
    int 		  lowest_table_dif_value;	/**< Lowest value in the diff matrix. */
    unsigned int word_size_k;				/**< Word size. */
    unsigned int fq_word_count;				/**< Total read words. */
    unsigned int ref_word_count;			/**< Total words of the reference sequence. */
    double norm;					/**< Value to normalize output matrix. */ 
    int dim_n;						/**< Dimension of the tables (dim_n X dim_n). */
    int memory_size;					/**< Memory size for each table. */
    double f_x;						/**< X coordinate of the coloured point. */
    double f_y;						/**< Y coordinate of the coloured point. */
    int base_quality;
    char pgm_fastq_filename[MAX_PGM_FILENAME_LENGTH];	/**< The three name of the files are stored for copy to the qc_report. */
    char pgm_quality_filename[MAX_PGM_FILENAME_LENGTH];
    char pgm_diff_filename[MAX_PGM_FILENAME_LENGTH];  
} chaos_game_data_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Initializes a Genomic Signature header
*  @param[in,out] header_gs_p pointer to the structure to fill
*  @param gs_filename genomic signature filename
*  @param k_cg word size
*  @param dim X and Y dimension
*  @return void
*  
*  This function initializes a Genomic Signature header with the header content
*  of the specified filename
*/
void header_gs_init(header_gs_t* header_gs_p, char* gs_filename, int k_cg, int dim);

/**
*  @brief Initializes the look-up table with nt codification
*  @param[in,out] lut_table pointer to the look-up table to fill
*  @return void
*  
*  This function initializes the look-up table with nt codification
*/
void lut_table_init(char* lut_table);

/**
*  @brief Initializes the genomic signature matrix 
*  @param table_gs_p pointer to the matrix to initialize
*  @param dim_n dimension of the matrix
*  @return genomic signature matrix
*  
*  This function initializes the genomic signature matrix, memory 
*  allocation is performed for the structure
*/
unsigned int** chaos_game_table_gs_init(unsigned int** table_gs_p, int dim_n);

/**
*  @brief Frees the genomic signature matrix 
*  @param table_gs_p pointer to the matrix to initialize
*  @param dim_n dimension of the matrix
*  @return void
*  
*  This function frees all positions of the genomic signature matrix
*/
void chaos_game_table_gs_free(unsigned int** table_gs_p, int dim_n);

/**
*  @brief Initializes chaos game data
*  @param[in,out] chaos_game_data_p pointer to the chaos game data structure
*  @param table_gs_p pointer to the genomic signature table
*  @param k_cg word size in the chaos game
*  @param dim_n dimension of the matrices
*  @param base_quality base quality for quality values normalization
*  @return void
*  
*  This function initializes chaos game data using word size, base quality,
*  dimension of the matrices and genomic signature matrix*  
*/
void chaos_game_data_init(chaos_game_data_t* chaos_game_data_p, unsigned int** table_gs_p, int k_cg, int dim_n, int base_quality);

/**
*  @brief Frees the chaos game data
*  @param chaos_game_data_p pointer to the chaos game data structure
*  @return void
*  
*  This function frees all positions of the genomic signature matrix
*/
void chaos_game_data_free(chaos_game_data_t* chaos_game_data_p);

/**
*  @brief Fill tables of chaos game data structure
*  @param[in,out] chaos_game_data_p pointer to the chaos game data structure
*  @param batch_p pointer to the status batch to read
*  @param mode mode that determines if only valid reads or all reads compute in the chaos game
*  @return void
*  
*  This function fills chaos game tables by processing fastq reads and 
*  computing the calculation of matrix coordinates
*/
void chaos_game_fill_tables(chaos_game_data_t* chaos_game_data_p, status_batch_t* batch_p, int mode);

/**
*  @brief Loads genomic signature matrix of the reference genome into the chaos game data
*  @param[in,out] chaos_game_data_p pointer to the chaos game data structure
*  @param gs_filename genomic signature filename
*  @return void
*  
*  This function loads genomic signature matrix of the reference genome into the chaos 
*  game data with the content of the supplied genomic signature file
*/
void chaos_game_load_table_gs(chaos_game_data_t* chaos_game_data_p, char* gs_filename);

/**
*  @brief Loads genomic signature matrix in a matrix structure
*  @param table_gs_p pointer to the matrix where data will be loaded
*  @param dim_n X and Y dimensions of the matrix
*  @param gs_filename genomic signature filename
*  @return genomic signature matrix
*  
*  This function loads genomic signature matrix of the reference genome into the chaos 
*  game data with the content of the supplied genomic signature file
*/
unsigned int chaos_game_load_table_gs_direct(unsigned int** table_gs_p, int dim_n, char* gs_filename);

/**
*  @brief Calculates the substraction between chaos game matrices
*  @param [in,out] chaos_game_data_p pointer to the chaos game data structure
*  @return void
*  
*  This function calculates the difference between chaos game matrix of the reference genome 
*  and the same matrix calculated from fastq file reads. The difference is performed by normalizing 
*  both matrices and substract values in the matrix
*/
void chaos_game_calculate_table_dif(chaos_game_data_t* chaos_game_data_p);

/**
*  @brief Validates the difference matrix
*  @param chaos_game_data_p pointer to the chaos game data structure
*  @return 1 if validated, 0 if not validated
*  
*  This function calculates mean value and standard deviation of the difference matrix
*  and validates both values. In the actual implementation all values are considered correct
*/
int chaos_game_validate_table_dif(chaos_game_data_t* chaos_game_data_p);

/**
*  @brief Prints to disk the content of the chaos game matrices
*  @param chaos_game_data_p pointer to the chaos game data structure
*  @param fq_path 
*  @param report_directory 
*  @return void
*  
*  This function writes matrices content into disk in the specified report directory.
*  Names of the file are made by adding suffixes to fastq filename
*/
void chaos_game_write_table_images(chaos_game_data_t* chaos_game_data_p, char* fq_path, char* report_directory);

/**
*  @brief Prints a matrix by console
*  @param table pointer to matrix
*  @param dim_n X and Y dimension of the matrix
*  @return void
*  
*  This function prints values of a square matrix with chaos game values by console for 
*  debug purposes
*/
void chaos_game_print_table(unsigned int** table, int dim_n);

#endif	/*  CHAOS_GAME_H  */



