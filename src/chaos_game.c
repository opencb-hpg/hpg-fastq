
#include <omp.h>

#include "chaos_game.h"

/* **********************************************
 *    		Private functions  		*
 * *********************************************/

void chaos_game_normalize_quality_table_(chaos_game_data_t* chaos_game_data_p);
void chaos_game_absolute_diff_table_(chaos_game_data_t* chaos_game_data_p);
void chaos_game_generate_pgm_file_(unsigned int** table, int dim_n, double norm, int k_cg, char* pgm_filename);


/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

/* **********************************************************************
 *    		Genomic signature header management functions  		*
 * *********************************************************************/

void header_gs_init(header_gs_t* header_gs_p, char* gs_filename, int k_cg, int dim) {
    strncpy(header_gs_p->gs_filename, gs_filename, 180);
    header_gs_p->word_size_k = k_cg;
    header_gs_p->dim_x = dim;
    header_gs_p->dim_x = dim;
    header_gs_p->ref_word_count = 0;
}

void lut_table_init(char* lut_table) {
    for (int i = 0; i < 256; i++) {
        lut_table[i] = NUM_LETTERS;
    }

    lut_table['A'] = 0;
    lut_table['a'] = 0;

    lut_table['C'] = 1;
    lut_table['c'] = 1;

    lut_table['G'] = 2;
    lut_table['g'] = 2;

    lut_table['T'] = 3;
    lut_table['t'] = 3;

    lut_table['N'] = N_UNDEFINED_BASE;
    lut_table['n'] = N_UNDEFINED_BASE;

    lut_table['>'] = GREATER_THAN;
}

/* **************************************************************
 *    		chaos_game_data management functions  		*
 * *************************************************************/

unsigned int** chaos_game_table_gs_init(unsigned int** table_gs_p, int dim_n) {
    table_gs_p = (unsigned int**) calloc(dim_n, sizeof(unsigned int*));

    int i;
    for (i = 0; i < dim_n; i++) {
        table_gs_p[i] = (unsigned int*) calloc(dim_n, sizeof(int));
    }

    return table_gs_p;
}

void chaos_game_table_gs_free(unsigned int** table_gs_p, int dim_n) {
    for (int i = 0; i < dim_n; i++) {
        free(table_gs_p[i]);
    }

    free(table_gs_p);
}

void chaos_game_data_init(chaos_game_data_t* chaos_game_data_p, unsigned int** table_gs_p, int k_cg, int dim_n, int base_quality) {
    int size_of_int, size_of_unsigned_int_pointer;

    chaos_game_data_p->word_size_k = k_cg;
    chaos_game_data_p->ref_word_count = 0; //Incremented later
    chaos_game_data_p->norm = 0;   //Calculated later with final count value

    chaos_game_data_p->dim_n = 1 << k_cg;    //dim_n = 2^k
    chaos_game_data_p->memory_size = 1 << (2 * k_cg);  //memory_size = 4^k = 2^(2*k)

    chaos_game_data_p->f_x = (double)(chaos_game_data_p->dim_n  * 0.5);
    chaos_game_data_p->f_y = chaos_game_data_p->f_x;
    chaos_game_data_p->base_quality = base_quality;

    chaos_game_data_p->mean_table_dif_value = 0;
    chaos_game_data_p->standard_deviation_table_dif_value = 0;
    chaos_game_data_p->highest_table_dif_value = 0;
    chaos_game_data_p->lowest_table_dif_value = 0;

    dim_n = chaos_game_data_p->dim_n;
    size_of_int = sizeof(int);
    size_of_unsigned_int_pointer = sizeof(unsigned int*);

    chaos_game_data_p->table_seq = (unsigned int**) calloc(dim_n, size_of_unsigned_int_pointer);
    chaos_game_data_p->table_q = (unsigned int**) calloc(dim_n, size_of_unsigned_int_pointer);
    chaos_game_data_p->table_gs = table_gs_p;
    chaos_game_data_p->table_dif = (unsigned int**) calloc(dim_n, size_of_unsigned_int_pointer);

    for (int i = 0; i < dim_n; i++) {
        chaos_game_data_p->table_seq[i] = (unsigned int*) calloc(dim_n, size_of_int);
        chaos_game_data_p->table_q[i] = (unsigned int*) calloc(dim_n, size_of_int);
        chaos_game_data_p->table_dif[i] = (unsigned int*) calloc(dim_n, size_of_int);
    }
}

void chaos_game_data_free(chaos_game_data_t* chaos_game_data_p) {
    if (chaos_game_data_p == NULL) return;

    int dim_n = chaos_game_data_p->dim_n;

    for (int i = 0; i < dim_n; i++) {
        free(chaos_game_data_p->table_seq[i]);
        free(chaos_game_data_p->table_q[i]);
        free(chaos_game_data_p->table_dif[i]);
    }

    if (chaos_game_data_p->table_seq != NULL) {
        free(chaos_game_data_p->table_seq);
        chaos_game_data_p->table_seq = NULL;
    }

    if (chaos_game_data_p->table_q != NULL) {
        free(chaos_game_data_p->table_q);
        chaos_game_data_p->table_q = NULL;
    }

    if (chaos_game_data_p->table_dif != NULL) {
        free(chaos_game_data_p->table_dif);
        chaos_game_data_p->table_dif = NULL;
    }

    free(chaos_game_data_p);
}

/* **********************************************************************
 *    		Genomic signature generation functions  		*
 * *********************************************************************/

void chaos_game_fill_tables(chaos_game_data_t* chaos_game_data_p, status_batch_t* batch_p, int mode) {
    int dim_n, word_size, base_quality, word_quality_substract, nt_word_count = 0;
    int num_reads, read_position, quality_position, read_length;
    unsigned int acc_word_quality = 0;
    int co_x, co_y;
    double f_x, f_y;
    char quality_character;
    char* seq;
    char* quality;
    fastq_batch_t* read_p;

    read_p = batch_p->read_p;
    seq = batch_p->read_p->seq;
    quality = batch_p->read_p->quality;
    num_reads = read_p->num_reads;
    f_x = chaos_game_data_p->f_x;
    f_y = chaos_game_data_p->f_y;
    dim_n = chaos_game_data_p->dim_n;
    word_size = chaos_game_data_p->word_size_k;
    base_quality = chaos_game_data_p->base_quality;
    word_quality_substract = base_quality * word_size;

    for (int i = 0; i < num_reads; i++) {
        if ((mode == ONLY_VALID_READS) && (batch_p->read_status[i] != VALID_READ)) {
            continue;
        }

        read_position = read_p->data_indices[i];
        read_length = read_p->data_indices[i+1] - read_p->data_indices[i];
        quality_position = read_position;

        for (int j = 0; j < read_length; j++) {

            quality_character = quality[quality_position++];

            switch (seq[read_position++]) {
                case 65: {   //"A"
                    f_x = f_x + ((dim_n - f_x) * 0.5);
                    f_y = f_y * 0.5;
                    nt_word_count++;
                    acc_word_quality += quality_character;
                    break;
                };
                case 67: {   //"C"
                    f_x = f_x * 0.5;
                    f_y = f_y * 0.5;
                    nt_word_count++;
                    acc_word_quality += quality_character;
                    break;
                };
                case 71: {   //"G"
                    f_x = f_x * 0.5;
                    f_y = f_y + ((dim_n - f_y) * 0.5);
                    nt_word_count++;
                    acc_word_quality += quality_character;
                    break;
                };
                case 84: {   //"T"
                    f_x = f_x + ((dim_n - f_x) * 0.5);
                    f_y = f_y + ((dim_n - f_y) * 0.5);
                    nt_word_count++;
                    acc_word_quality += quality_character;
                    break;
                };
                case 78: {
                    nt_word_count = 0;
                    acc_word_quality = 0;
                    break;
                }
            }

            if (nt_word_count == word_size) { //a whole word is completed, compute coordinate values
                co_x = (int) f_x;
                co_y = (int) f_y;


                //if coordinate is on the boundary the integer coordinate is decremented 1
                //and the double coordinate is decremented an epsilon infinitesimal value
                if (co_x == dim_n)  {
                    co_x = dim_n - 1;
                    f_x = f_x - EPSILON;
                }

                if (co_y == dim_n)  {
                    co_y = dim_n - 1;
                    f_y = f_y - EPSILON;
                }

                chaos_game_data_p->table_seq[co_x][co_y]++;  //only when nts in work complete a whole word
                chaos_game_data_p->fq_word_count++;
                nt_word_count--;

                chaos_game_data_p->table_q[co_x][co_y] += (acc_word_quality - word_quality_substract); //Substract base_quality * word_size

                acc_word_quality = acc_word_quality - quality[quality_position - word_size];
            }
        } //end of for loop, read finished

        nt_word_count = 0;
        acc_word_quality = 0;
    } //end of for loop, batch finished

}

void chaos_game_load_table_gs(chaos_game_data_t* chaos_game_data_p, char* gs_filename) {
    FILE* gs_fd;  //genomic signature file descriptor
    header_gs_t* header_gs = (header_gs_t*) calloc(1, sizeof(header_gs_t));

    int dim_n = chaos_game_data_p->dim_n;

    if ((gs_fd = fopen(gs_filename, "rb")) == NULL) {
        LOG_FATAL("Cannot open Genomic Signature file for read, aborting execution\n");
    }

    fread(header_gs, 1, sizeof(header_gs_t), gs_fd);
    strncpy(header_gs->gs_filename, gs_filename, 180);
    chaos_game_data_p->ref_word_count = header_gs->ref_word_count;

// printf("gs_filename: %s ............\n", gs_filename);
// printf("chaos_game_data_p->fq_word_count: %i\n", chaos_game_data_p->fq_word_count);
// printf("header_gs->dim_x: %u, header_gs->dim_y: %u\n", header_gs->dim_x, header_gs->dim_y);
// printf("header_gs->gs_filename: %s\n", header_gs->gs_filename);
// printf("header_gs->ref_word_count: %u\n", header_gs->ref_word_count);
// printf("header_gs->word_size_k: %u\n", header_gs->word_size_k);

    for (int i = 0; i < dim_n; i++) {
        fread(chaos_game_data_p->table_gs[i], dim_n, sizeof(int), gs_fd);
    }

    free(header_gs);
    fclose(gs_fd);
}

unsigned int chaos_game_load_table_gs_direct(unsigned int** table_gs_p, int dim_n, char* gs_filename) {
    FILE* gs_fd;  //genomic signature file descriptor
    int ref_word_count;
    header_gs_t* header_gs = (header_gs_t*) calloc(1, sizeof(header_gs_t));

    if ((gs_fd = fopen(gs_filename, "rb")) == NULL) {
        LOG_FATAL("Cannot open Genomic Signature file for read, aborting execution\n");
    }

    fread(header_gs, 1, sizeof(header_gs_t), gs_fd);

    strncpy(header_gs->gs_filename, gs_filename, 180);
    ref_word_count = header_gs->ref_word_count;

    for (int i = 0; i < dim_n; i++) {
        fread(table_gs_p[i], dim_n, sizeof(int), gs_fd);
    }

    free(header_gs);
    fclose(gs_fd);

    return ref_word_count;
}

void chaos_game_calculate_table_dif(chaos_game_data_t* chaos_game_data_p) {
    int highest_value, lowest_value;
    double fq_norm, gs_norm;  //fastq norm and genomic signature norm
    char log_message[100];

    int dim_n = chaos_game_data_p->dim_n;
    int memory_size = chaos_game_data_p->memory_size;

    //calculation of fastq norm value
    fq_norm = (double)(1.0 * chaos_game_data_p->fq_word_count / memory_size);

    if (fq_norm > 0.0) {
        fq_norm = (128.0 / fq_norm);
    } else {
        LOG_FATAL("Error in fastq norm calculation, aborting program\n");
    }

    sprintf(log_message, "Fastq file normalization ratio = %12.6f\n", fq_norm);
    LOG_DEBUG(log_message);

    //calculation of genomic signature value
    gs_norm = (double)(1.0 * chaos_game_data_p->ref_word_count / memory_size);

    if (gs_norm > 0.0) {
        gs_norm = (128.0 / gs_norm);
    } else {
        LOG_FATAL("Error in genomic signature norm calculation, aborting program\n");
    }

    sprintf(log_message, "Genomic signature normalization ratio = %12.6f\n", gs_norm);
    LOG_DEBUG(log_message);

    //initial values
    highest_value = -32768;
    lowest_value = 32768;

    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {
            chaos_game_data_p->table_dif[i][j] = (chaos_game_data_p->table_seq[i][j] * fq_norm) - (chaos_game_data_p->table_gs[i][j] * gs_norm);
	    
            if (highest_value < chaos_game_data_p->table_dif[i][j]) highest_value = chaos_game_data_p->table_dif[i][j];
            if (lowest_value > chaos_game_data_p->table_dif[i][j]) lowest_value = chaos_game_data_p->table_dif[i][j];
        }
    }

    chaos_game_data_p->highest_table_dif_value = highest_value;
    chaos_game_data_p->lowest_table_dif_value = lowest_value;

    sprintf(log_message, "Interval of variation of diff matrix values = [%d, %d]\n", highest_value, lowest_value);
    LOG_DEBUG(log_message);
}

int chaos_game_validate_table_dif(chaos_game_data_t* chaos_game_data_p) {
    int dim_n = chaos_game_data_p->dim_n;

    double mean_table_dif_value;
    double standard_deviation_table_dif_value;

    //calculation of mean value with positive and negative value
    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {
            mean_table_dif_value += chaos_game_data_p->table_dif[i][j];
        }
    }

    mean_table_dif_value /= (double)(1.0 * dim_n * dim_n);

    //calculation of the standard deviation    
    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {
            standard_deviation_table_dif_value += pow((chaos_game_data_p->table_dif[i][j] - mean_table_dif_value), 2);
        }
    }

    standard_deviation_table_dif_value /= (double)(1.0 * dim_n * dim_n);
    standard_deviation_table_dif_value = sqrt(standard_deviation_table_dif_value);

    chaos_game_data_p->mean_table_dif_value = mean_table_dif_value;
    chaos_game_data_p->standard_deviation_table_dif_value = standard_deviation_table_dif_value;

    //TODO: return 0 or 1 using heuristic values with known precaculated thresholds

    return 1;
}

void chaos_game_write_table_images(chaos_game_data_t* chaos_game_data_p, char* fq_path, char* report_directory) {
    char pgm_base_filename[MAX_PGM_FILENAME_LENGTH];
    char pgm_filename[MAX_PGM_FILENAME_LENGTH];
    char fq_filename[MAX_PGM_FILENAME_LENGTH];
    char log_message[200];
    double fq_norm, q_norm;  //Normalization factors for fastq sequences, qualities and genomic signature of the reference sequence

    int dim_n = chaos_game_data_p->dim_n;

    get_filename_from_path(fq_path, fq_filename);
    sprintf(pgm_base_filename, "%s/%s%s%d", report_directory, fq_filename, K_VALUE_INFIX, chaos_game_data_p->word_size_k);

    //Genomic Signature for fastq sequences
    sprintf(pgm_filename, "%s%s", pgm_base_filename, FASTQ_PGM_FILENAME_SUFFIX);
    strcpy(chaos_game_data_p->pgm_fastq_filename, pgm_filename);

    sprintf(log_message, "Words read in FastQ_Source_File = %s: %u\n", fq_filename, chaos_game_data_p->fq_word_count);
    LOG_DEBUG(log_message);
    sprintf(log_message, "Genomic Signature of fastq in file: %s\n", pgm_filename);
    LOG_DEBUG(log_message);

    fq_norm = (double) chaos_game_data_p->fq_word_count;
    fq_norm = fq_norm / chaos_game_data_p->memory_size;

    if (fq_norm > 0.0) {
        fq_norm = 128.0 / fq_norm;
    } else {
        LOG_FATAL("Error in fastq sequences normalization, aborting execution\n");
    }

    chaos_game_generate_pgm_file_(chaos_game_data_p->table_seq, dim_n, fq_norm, chaos_game_data_p->word_size_k, pgm_filename);

    //Genomic Signature for qualities
    sprintf(pgm_filename, "%s%s", pgm_base_filename, QUALITY_PGM_FILENAME_SUFFIX);
    strcpy(chaos_game_data_p->pgm_quality_filename, pgm_filename);

    sprintf(log_message, "Genomic Signature of qualities in file: %s\n", pgm_filename);
    LOG_DEBUG(log_message);

    //Quality values for each word must be normalized to its respective word frequency and to the work size
    chaos_game_normalize_quality_table_(chaos_game_data_p);

    q_norm = 256.0 / MAX_QUALITY_IN_TABLE;

    chaos_game_generate_pgm_file_(chaos_game_data_p->table_q, dim_n , q_norm, chaos_game_data_p->word_size_k, pgm_filename);


    //Genomic Signature for difference matrix between fastq sequences and reference genome
    sprintf(pgm_filename, "%s%s", pgm_base_filename, DIFF_PGM_FILENAME_SUFFIX);
    strcpy(chaos_game_data_p->pgm_diff_filename, pgm_filename);

    sprintf(log_message, "Genomic Signature of diff in file: %s\n", pgm_filename);
    LOG_DEBUG(log_message);

    chaos_game_absolute_diff_table_(chaos_game_data_p);

    //cast to unsigned int without no conflict, all values are greater than 0
    chaos_game_generate_pgm_file_((unsigned int**) chaos_game_data_p->table_dif, dim_n, 1, chaos_game_data_p->word_size_k, pgm_filename);
}

void chaos_game_print_table(unsigned int** table, int dim_n) {
    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {
            printf("table[%i][%i]: %u\n", i, j, table[i][j]);
        }
    }
}

/* **************************************************************
 *    		Private functions implementations  		*
 * *************************************************************/

void chaos_game_normalize_quality_table_(chaos_game_data_t* chaos_game_data_p) {
    int dim_n = chaos_game_data_p->dim_n;
    int word_size = chaos_game_data_p->word_size_k;

    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {
            if (chaos_game_data_p->table_seq[i][j] > 0) {  //Avoid division by 0
                chaos_game_data_p->table_q[i][j] /= word_size;     //Normalization to word_size
                chaos_game_data_p->table_q[i][j] /= chaos_game_data_p->table_seq[i][j]; //Normalization to word frequency
            } else {
                chaos_game_data_p->table_q[i][j] = 0.0;     //0 if no word detected
            }
        }
    }
}

void chaos_game_absolute_diff_table_(chaos_game_data_t* chaos_game_data_p) {
    int dif_value;
    int dim_n = chaos_game_data_p->dim_n;

    for (int i = 0; i < dim_n; i++) {
        for (int j = 0; j < dim_n; j++) {

            dif_value = abs(chaos_game_data_p->table_dif[i][j]);

            if (dif_value > 255) {
                dif_value = 255;
            }

            chaos_game_data_p->table_dif[i][j] = (uchar) dif_value; //normalized to expected values
        }
    }
}

void chaos_game_generate_pgm_file_(unsigned int** table, int dim_n, double norm, int k_cg, char* pgm_filename) {
    int redim_n, int_point_value, zoom_factor, dif_k;
    float float_point_value;
    uchar** image_matrix;  //for output images
    FILE* pgm_fd;

    //if k for chaos game is less than MIN_K_IMAGE_VALUE image must be resized (zoom)
    if (k_cg < MIN_K_IMAGE_VALUE) {
        redim_n = MIN_IMAGE_PIXEL_SIZE;  //Default minimum: 128 X 128
    } else {
        redim_n = dim_n;
    }

    //memory allocation for image matrix
    image_matrix =  calloc(redim_n, sizeof(uchar *));    //initialized to zero by calloc

    for (int i = 0; i < redim_n; i++) {
        image_matrix[i] =  calloc(redim_n, sizeof(uchar));
    }

    //.pgm file is opened for write in BINARY mode
    pgm_fd = fopen(pgm_filename, "wb");

    if (k_cg >= MIN_K_IMAGE_VALUE) {
        for (int i = 0; i < dim_n; i++) {
            for (int j = 0; j < dim_n; j++) {
                float_point_value = (float) table[i][j];
                float_point_value =  float_point_value * norm;
                int_point_value = (int) float_point_value;
                image_matrix[i][j] = (uchar) int_point_value;
            }
        }
    } else {
        dif_k = MIN_K_IMAGE_VALUE - k_cg;
        zoom_factor = pow(2, dif_k); //zoom_factor = 2 ^ dif_k
        
        //zoom_factor = 2 ^ dif_k
//         zoom_factor = 1;
//         for (int i = 0; i < dif_k; i++) {
//             zoom_factor *= 2; 
//         }

        for (int i = 0; i < dim_n; i++) {
            for (int j = 0; j < dim_n; j++) {
                float_point_value = (float) table[i][j];
                float_point_value =  float_point_value * norm;
                int_point_value  = (int) float_point_value;

                for (int ii = 0; ii < zoom_factor; ii++) {
                    for (int jj = 0; jj < zoom_factor; jj++) {
                        image_matrix[(i * zoom_factor) + ii][(j * zoom_factor)+jj] = (uchar) int_point_value;
                    }
                }

            }
        }
    }

    //write image data to file
    fprintf(pgm_fd, "P5\n");
    fprintf(pgm_fd, "%d %d\n", redim_n, redim_n); //matrix dimensions
    fprintf(pgm_fd, "%d\n", 255);        //number of grey tones

    for (int i = 0; i < redim_n; i++) {   //redim_n >= dim_n with every condition
        fwrite(image_matrix[i], redim_n, 1, pgm_fd);
    }

    printf("\n");

    if (k_cg == MIN_K_IMAGE_VALUE) {
        for (int i = 0; i < dim_n; i++) {
            free(image_matrix[i]);
        }
    } else {
        for (int i = 0; i < redim_n; i++) {
            free(image_matrix[i]);
        }
    }

    free(image_matrix);
}