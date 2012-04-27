
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "qc_report.h"

/* **********************************************
 *    		Private functions  		*
 * *********************************************/

void print_qc_text_report_file(qc_report_t qc_report, int base_quality, char *inputfilename, char *outfilename);
void print_chaos_game_text_report_file(qc_report_t qc_report, char *outfilename);
void generate_position_datafile(qc_report_t qc_report, char* filename);
void generate_read_quality_datafile(qc_report_t* qc_report_p, char* quality_filename);
void generate_kmers_datafile(qc_report_t* qc_report_p, char* kmers_filename);
void generate_kmers_position_datafile(qc_report_t* qc_report_p, char* kmers_position_filename);
qc_graph_t* set_qc_graph_default_values(qc_graph_t* qc_graph);
void plot_nt_quality(qc_report_t* qc_report_p, char* datafilename, char* graphfilename);
void plot_read_quality(char* quality_filename, char* graph_filename);
void plot_sequence_length_distribution(char* data_filename, char* graph_filename);
void plot_base_sequence_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
void plot_base_gc_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
void generate_gc_histogram_datafile(qc_report_t* qc_report_p, char* gc_histogram_filename);
void plot_sequence_gc_content(char* data_filename, char* graph_filename);
void plot_n_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
void plot_kmers_position(qc_report_t* qc_report_p, char* data_filename, char* graph_filename);
void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename);
void generate_html_file_with_images(qc_report_t qc_report, int kmers_on, char* html_filename, char* report_directory, char* in_shortname);
void generate_html_file_valid_with_images(qc_report_t qc_report, int kmers_on, int cg_flag, char* html_filename, char* report_directory, char* in_shortname, int valid);

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

void qc_report_get_values_from_chaos_game(chaos_game_data_t* chaos_game_data_p, qc_report_t* qc_report) {
    strcpy(qc_report->pgm_fastq_filename, chaos_game_data_p->pgm_fastq_filename);
    strcpy(qc_report->pgm_quality_filename, chaos_game_data_p->pgm_quality_filename);
    strcpy(qc_report->pgm_diff_filename, chaos_game_data_p->pgm_diff_filename);
    qc_report->dim_n  = chaos_game_data_p->dim_n;
    qc_report->word_size_k = chaos_game_data_p->word_size_k;
    qc_report->fq_word_count = chaos_game_data_p->fq_word_count;
    qc_report->mean_table_dif_value = chaos_game_data_p->mean_table_dif_value;
    qc_report->standard_deviation_table_dif_value = chaos_game_data_p->standard_deviation_table_dif_value;
    qc_report->highest_table_dif_value = chaos_game_data_p->highest_table_dif_value;
    qc_report->lowest_table_dif_value = chaos_game_data_p->lowest_table_dif_value;
}

void generate_report(qc_report_t qc_report, char* inputfilename, int base_quality, int kmers_on, int cg_flag, char* report_directory, int valid) {  
    //if there is no read there is no report to print!!!!
    if (qc_report.nb_reads == 0) return;

    char in_shortname[MAX_FULL_PATH_LENGTH];
    char processed_filename[MAX_FULL_PATH_LENGTH];
    char data_filename[MAX_FULL_PATH_LENGTH];
    char quality_filename[MAX_FULL_PATH_LENGTH];
    char gc_histogram_filename[MAX_FULL_PATH_LENGTH];
    char graph_filename[MAX_FULL_PATH_LENGTH];
    char kmers_filename[MAX_FULL_PATH_LENGTH];
    char kmers_position_filename[MAX_FULL_PATH_LENGTH];
    char html_filename[MAX_FULL_PATH_LENGTH];
    char* str_valid_suffix = "";

    get_filename_from_path(inputfilename, in_shortname);

    str_valid_suffix = (valid) ? (char*) ".valid" : (char*)  ".invalid";

    // print qc text report
    sprintf(data_filename, "%s/%s%s%s", report_directory, in_shortname, QC_SUFFIX, str_valid_suffix);
    sprintf(processed_filename, "%s%s", in_shortname, str_valid_suffix);

    print_qc_text_report_file(qc_report, base_quality, processed_filename, data_filename);

    // print chaos game text report
    if (cg_flag) {
        sprintf(data_filename, "%s/%s%s%s", report_directory, in_shortname, CG_SUFFIX, str_valid_suffix);
        print_chaos_game_text_report_file(qc_report, data_filename);
    }

    // generate datafile and plots for postions
    sprintf(data_filename, "%s/%s%s%s", report_directory, in_shortname, POSITION_DATA_SUFFIX, str_valid_suffix);
    generate_position_datafile(qc_report, data_filename);

    sprintf(graph_filename, "%s/%s.nt_quality%s.png", report_directory, in_shortname, str_valid_suffix);
    plot_nt_quality(&qc_report, data_filename, graph_filename);

    sprintf(quality_filename, "%s/%s%s%s", report_directory, in_shortname, QUALITY_DATA_SUFFIX, str_valid_suffix);
    sprintf(graph_filename, "%s/%s.quality_read%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_read_quality_datafile(&qc_report, quality_filename);
    plot_read_quality(quality_filename, graph_filename);

    sprintf(graph_filename, "%s/%s.read_length%s.png", report_directory, in_shortname, str_valid_suffix);
    plot_sequence_length_distribution(data_filename, graph_filename);

    sprintf(graph_filename, "%s/%s.base_sequence%s.png", report_directory, in_shortname, str_valid_suffix);
    plot_base_sequence_content(&qc_report, data_filename, graph_filename);

    sprintf(graph_filename, "%s/%s.base_gc%s.png", report_directory, in_shortname, str_valid_suffix);
    plot_base_gc_content(&qc_report, data_filename, graph_filename);

    sprintf(gc_histogram_filename, "%s/%s%s%s", report_directory, in_shortname, GC_DATA_SUFFIX, str_valid_suffix);
    sprintf(graph_filename, "%s/%s.gc_histogram%s.png", report_directory, in_shortname, str_valid_suffix);
    generate_gc_histogram_datafile(&qc_report, gc_histogram_filename);
    plot_sequence_gc_content(gc_histogram_filename, graph_filename);

    if (kmers_on) {
        sprintf(kmers_filename, "%s/%s%s%s", report_directory, in_shortname, KMERS_TOTAL_SUFFIX, str_valid_suffix);
        sprintf(kmers_position_filename, "%s/%s%s%s", report_directory, in_shortname, KMERS_POSITION_SUFFIX, str_valid_suffix);
        sprintf(graph_filename, "%s/%s.kmers_positions%s.png", report_directory, in_shortname, str_valid_suffix);
        generate_kmers_datafile(&qc_report, kmers_filename);
        generate_kmers_position_datafile(&qc_report, kmers_position_filename);
        plot_kmers_position(&qc_report, kmers_position_filename, graph_filename);
    }

    sprintf(graph_filename, "%s/%s.n_base%s.png", report_directory, in_shortname, str_valid_suffix);
    plot_n_content(&qc_report, data_filename, graph_filename);

    sprintf(html_filename, "%s/%s%s%s", report_directory, in_shortname, str_valid_suffix, HTML_FILE_SUFFIX);
    generate_html_file_valid_with_images(qc_report, kmers_on, cg_flag, html_filename, report_directory, in_shortname, valid);
}

void print_qc_text_report_file(qc_report_t qc_report, int base_quality, char *inputfilename, char *outfilename) {
    FILE* fd = (outfilename == NULL ? stdout : fopen(outfilename, "w"));

    fprintf(fd, "\n----------------------------------------------\n");
    fprintf(fd, "            Q C      R E P O R T                  ");
    fprintf(fd, "\n----------------------------------------------\n");
    fprintf(fd, "\nProcessed file  : %s\n", inputfilename);
    fprintf(fd, "\nNumber of reads  : %li\n", qc_report.nb_reads);
    fprintf(fd, "\nMin. read length : %i\n", qc_report.min_read_length);
    fprintf(fd, "\nMax. read length : %i\n", qc_report.max_read_length);
    fprintf(fd, "\nMean read length : %i\n", qc_report.mean_read_length);
    fprintf(fd, "\nMean read quality: %li [%c]\n", qc_report.mean_read_quality, (qc_report.mean_read_quality + base_quality));

    fprintf(fd, "\nAnaylisis of nucleotides (A, C, T, G, N):\n");
    fprintf(fd, "\tA: %02.2f %\n", qc_report.a_perc);
    fprintf(fd, "\tC: %02.2f %\n", qc_report.c_perc);
    fprintf(fd, "\tG: %02.2f %\n", qc_report.g_perc);
    fprintf(fd, "\tT: %02.2f %\n", qc_report.t_perc);
    fprintf(fd, "\tN: %02.2f %\n", qc_report.n_perc);

    fprintf(fd, "\nGC content:\n");
    fprintf(fd, "\tGC: %02.2f %\n", qc_report.c_perc + qc_report.g_perc);

    fprintf(fd, "\nMean quality per nucleotide position:\n");

    int i;
    for (i = 0; i < MAX_LINE_LENGTH; i = i + 5) {
        if (qc_report.mean_nt_quality[i] > 0) {
            fprintf(fd, "\tpos. %3i: %3li [%c]\t", i + 1, qc_report.mean_nt_quality[i], (qc_report.mean_nt_quality[i] + base_quality));
            fprintf(fd, "\tpos. %3i: %3li [%c]\t", i + 2, qc_report.mean_nt_quality[i+1], (qc_report.mean_nt_quality[i+1] + base_quality));
            fprintf(fd, "\tpos. %3i: %3li [%c]\t", i + 3, qc_report.mean_nt_quality[i+2], (qc_report.mean_nt_quality[i+2] + base_quality));
            fprintf(fd, "\tpos. %3i: %3li [%c]\t", i + 4, qc_report.mean_nt_quality[i+3], (qc_report.mean_nt_quality[i+3] + base_quality));
            fprintf(fd, "\tpos. %3i: %3li [%c]\n", i + 5, qc_report.mean_nt_quality[i+4], (qc_report.mean_nt_quality[i+4] + base_quality));
        }
    }

    if (outfilename != NULL)  fclose(fd);
}

void print_chaos_game_text_report_file(qc_report_t qc_report, char *outfilename) {
    FILE* fd = (outfilename == NULL ? stdout : fopen(outfilename, "w"));

    fprintf(fd, "\n---------------------------------------------------------------------\n");
    fprintf(fd, "       G E N O M I C   S I G N A T U R E   ( C H A O S   G A M E )       ");
    fprintf(fd, "\n---------------------------------------------------------------------\n");
    fprintf(fd, "\nFastq sequences image filename : %s\n", qc_report.pgm_fastq_filename);
    fprintf(fd, "\nFastq qualities image filename : %s\n", qc_report.pgm_quality_filename);
    fprintf(fd, "\nDiff image filename : %s\n", qc_report.pgm_diff_filename);
    fprintf(fd, "\nWord size (k) : %i\n", qc_report.word_size_k);
    fprintf(fd, "\nNumber of different words (k exp 4) : %i (%i X %i)\n", (qc_report.dim_n * qc_report.dim_n), qc_report.dim_n, qc_report.dim_n);
    fprintf(fd, "\nWords located in fastq sequences : %u\n", qc_report.fq_word_count);
    fprintf(fd, "\nMean value of the diff matrix : %.3f\n", qc_report.mean_table_dif_value);
    fprintf(fd, "\nStandard deviation of the diff matrix : %.3f\n", qc_report.standard_deviation_table_dif_value);
    fprintf(fd, "\nRange of variation of the diff matrix : [%i, %i]\n", qc_report.lowest_table_dif_value, qc_report.highest_table_dif_value);

    if (outfilename != NULL)  fclose(fd);
}

void generate_position_datafile(qc_report_t qc_report, char* filename) {
    FILE* fd = (filename == NULL ? stdout : fopen(filename, "w"));

    int nt_last_position = qc_report.nt_counter[0];
    int length, i;

    for (i = 1; i <= qc_report.max_read_length; i++) {

        if (nt_last_position != qc_report.nt_counter[i]) {
            length = nt_last_position - qc_report.nt_counter[i];
            nt_last_position = qc_report.nt_counter[i];
        } else {
            length = 0;
        }

        fprintf(fd, "%i\t%li\t%i\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n", i, qc_report.mean_nt_quality[i-1], length, (double)(100.0 * qc_report.nt_type_counter[i-1][A]) / qc_report.nt_counter[i-1], (double)(100.0 * qc_report.nt_type_counter[i-1][C]) / qc_report.nt_counter[i-1], (double)(100.0 * qc_report.nt_type_counter[i-1][G]) / qc_report.nt_counter[i-1], (double)(100.0 * qc_report.nt_type_counter[i-1][T]) / qc_report.nt_counter[i-1], (double)(100.0 * qc_report.nt_type_counter[i-1][N]) / qc_report.nt_counter[i-1], ((100.0 * (double) qc_report.nt_type_counter[i-1][C]) / qc_report.nt_counter[i-1] + (100.0 * (double) qc_report.nt_type_counter[i-1][G]) / qc_report.nt_counter[i-1]));
    }

    if (fd != NULL)  fclose(fd);
}

void generate_kmers_datafile(qc_report_t* qc_report_p, char* filename) {
    FILE* fd = (filename == NULL ? stdout : fopen(filename, "w"));

    fprintf(fd, "\n----------------------------------------------\n");
    fprintf(fd, "            K M E R S     C O U N T               ");
    fprintf(fd, "\n----------------------------------------------\n\n");

    int i;
    for (i = 0; i < KMERS_COMBINATIONS; i++) {
        //printf("qc_report_p->qc_kmers[%i].kmer: %s\n", i, qc_report_p->qc_kmers[i].kmer);
        fprintf(fd, "%s\t%i\n", qc_report_p->qc_kmers[i].kmer, qc_report_p->qc_kmers[i].total_count);
    }

    if (fd != NULL)  fclose(fd);
}

void generate_kmers_position_datafile(qc_report_t* qc_report_p, char* filename) {
    FILE* fd = (filename == NULL ? stdout : fopen(filename, "w"));

    int i, j;
    for (i = 0; i <= qc_report_p->max_read_length - 5; i++) {

        fprintf(fd, "%i", i);

        for (j = 0; j < KMERS_COMBINATIONS; j++) {
            fprintf(fd, "\t%i", qc_report_p->qc_kmers[j].position_count[i]);
        }

        fprintf(fd, "\n");
    }

    if (fd != NULL)  fclose(fd);
}

void plot_nt_quality(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Quality per nucleotide position";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Quality (normalized Prhed scale)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = NT_QUALITY_COLUMN;

printf("before generate_gnuplot_image - --------------->\n");
    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
printf("after generate_gnuplot_image - --------------->\n");    
}

//plot_nt_quality_public("/tmp/200k_reads_pe_1.fastq.pos.dat.valid", "/tmp/200k_reads_pe_1.fastq.nt_quality.valid.png");
// void plot_nt_quality_public(char* data_filename, char* graph_filename) {
//     qc_graph_t qc_graph;
//     set_qc_graph_default_values(&qc_graph);
// 
//     qc_graph.title = "Quality per nucleotide position";
//     qc_graph.xlabel = "nt position";
//     qc_graph.ylabel = "Quality (normalized Prhed scale)";
//     qc_graph.type = "lines";
//     qc_graph.x_autoscale = 0;
//     qc_graph.x_start = 0;
//     qc_graph.x_end = 100;
//     qc_graph.y_autoscale = 1;
//     qc_graph.x_column = POS_COLUMN;
//     qc_graph.num_y_columns = 1;
//     qc_graph.y_columns[0] = NT_QUALITY_COLUMN;
// 
// printf("before generate_gnuplot_image - --------------->\n");
//     generate_gnuplot_image(qc_graph, data_filename, graph_filename);
// printf("after generate_gnuplot_image - --------------->\n");    
// }

void plot_read_quality(char* quality_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per sequence quality scores";
    qc_graph.xlabel = "Quality (normalized Prhed scale)";
    qc_graph.ylabel = "Num. of reads";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = MAX_QUALITY_VALUE;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = READ_QUALITY_COLUMN;

    generate_gnuplot_image(qc_graph, quality_filename, graph_filename);
}

void generate_read_quality_datafile(qc_report_t* qc_report_p, char* quality_filename) {
    FILE* fd = (quality_filename == NULL ? stdout : fopen(quality_filename, "w"));

    int i;
    for (i = 0; i < MAX_QUALITY_VALUE; i++) {
        fprintf(fd, "%i\t%i\n", i, qc_report_p->per_sequence_quality[i]);
    }

    if (quality_filename != NULL)  fclose(fd);
}

void plot_sequence_length_distribution(char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Sequence length distribution";
    qc_graph.xlabel = "Sequence length";
    qc_graph.ylabel = "Num. of reads";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 1;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = NUM_READ_LENGTH_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
}

void plot_base_sequence_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base sequence content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Percentage of nt (%)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 4;
    qc_graph.y_columns[0] = A_COLUMN;
    qc_graph.y_columns[1] = C_COLUMN;
    qc_graph.y_columns[2] = G_COLUMN;
    qc_graph.y_columns[3] = T_COLUMN;
    qc_graph.y_titles[0] = "A %";
    qc_graph.y_titles[1] = "C %";
    qc_graph.y_titles[2] = "G %";
    qc_graph.y_titles[3] = "T %";

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
}

void plot_base_gc_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base GC content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "GC %";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = GC_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
}

void generate_gc_histogram_datafile(qc_report_t* qc_report_p, char* gc_histogram_filename) {
    FILE* fd = (gc_histogram_filename == NULL ? stdout : fopen(gc_histogram_filename, "w"));

    for (int i = 0; i <= 100; i++) {
        fprintf(fd, "%i\t%i\n", i, qc_report_p->gc_histogram[i]);
    }

    if (gc_histogram_filename != NULL)  fclose(fd);
}

void plot_sequence_gc_content(char* gc_histogram_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per sequence GC content";
    qc_graph.xlabel = "GC %";
    qc_graph.ylabel = "Num. of reads";
    qc_graph.type = "boxes";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = 99;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = GC_HISTOGRAM_COLUMN;

    generate_gnuplot_image(qc_graph, gc_histogram_filename, graph_filename);
}

void plot_n_content(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Per base sequence content";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Percentage of N (%)";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = POS_COLUMN;
    qc_graph.num_y_columns = 1;
    qc_graph.y_columns[0] = N_COLUMN;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
}

void plot_kmers_position(qc_report_t* qc_report_p, char* data_filename, char* graph_filename) {
    qc_graph_t qc_graph;
    set_qc_graph_default_values(&qc_graph);

    qc_graph.title = "Relative enrichment over read length";
    qc_graph.xlabel = "nt position";
    qc_graph.ylabel = "Num. of kmers";
    qc_graph.type = "lines";
    qc_graph.x_autoscale = 0;
    qc_graph.x_start = 0;
    qc_graph.x_end = qc_report_p->max_read_length;
    qc_graph.y_autoscale = 1;
    qc_graph.x_column = KMER_POS_COLUMN;
    qc_graph.num_y_columns = 5;
    qc_graph.y_columns[0] = KMER_1_COLUMN;
    qc_graph.y_columns[1] = KMER_2_COLUMN;
    qc_graph.y_columns[2] = KMER_3_COLUMN;
    qc_graph.y_columns[3] = KMER_4_COLUMN;
    qc_graph.y_columns[4] = KMER_5_COLUMN;
    qc_graph.y_titles[0] = qc_report_p->qc_kmers[0].kmer;
    qc_graph.y_titles[1] = qc_report_p->qc_kmers[1].kmer;
    qc_graph.y_titles[2] = qc_report_p->qc_kmers[2].kmer;
    qc_graph.y_titles[3] = qc_report_p->qc_kmers[3].kmer;
    qc_graph.y_titles[4] = qc_report_p->qc_kmers[4].kmer;

    generate_gnuplot_image(qc_graph, data_filename, graph_filename);
}

void get_reads_per_length_from_nt_counter(int nt_counter[], int* length_p) {
    int nt_last_position = nt_counter[0];

    for (int j = 0; j < MAX_LINE_LENGTH ; j++) {
        if (nt_last_position != nt_counter[j]) {
            length_p[j] = nt_last_position - nt_counter[j];
            nt_last_position = nt_counter[j];
        } else {
            length_p[j] = 0;
        }
    }
}

qc_graph_t* set_qc_graph_default_values(qc_graph_t* qc_graph) {
    qc_graph->x_autoscale = 1;
    qc_graph->x_start = 1;
    qc_graph->x_end = 100;
    qc_graph->y_autoscale = 1;
    qc_graph->y_start = 0;
    qc_graph->y_end = 100;
    qc_graph->lmargin = 10;
    qc_graph->rmargin = 4;
    qc_graph->tmargin = 3;
    qc_graph->bmargin = 4;
    qc_graph->title = "Title";
    qc_graph->xlabel = "X axis";
    qc_graph->ylabel = "Y axis";
    qc_graph->type = "lines";
    qc_graph->x_column = 0;
    qc_graph->num_y_columns = 1;
    qc_graph->y_columns[0] = 1;
    qc_graph->y_titles[0] = "";
    qc_graph->y_titles[1] = "";
    qc_graph->y_titles[2] = "";
    qc_graph->y_titles[3] = "";
    qc_graph->y_titles[4] = "";
    qc_graph->y_titles[5] = "";
    qc_graph->y_titles[6] = "";
    qc_graph->y_titles[7] = "";
    qc_graph->y_titles[8] = "";
    qc_graph->y_titles[9] = "";

    return qc_graph;
}

void generate_gnuplot_image(qc_graph_t qc_graph, char* data_filename, char* graph_filename) {
  
  //printf("data_filename: %s\n", data_filename);
  //printf("graph_filename: %s\n", graph_filename);

//   for (int i = 0; i < 10; i++) {
//       FILE* test_fd = popen("ls", "w");
//       //FILE* test_fd = fopen("/tmp/prueba.txt", "w");
//       printf("test_fd number %i open OK...\n", i);
//       pclose(test_fd);
//       //fclose(test_fd);
//       printf("test_fd number %i closed OK...\n", i);
//   }

  //FILE *graph_fd = popen("display", "w");
  //FILE *graph_fd = popen("gnuplot -persist", "w");
  //printf("closing graph_fd %x...\n", graph_fd);
  //pclose(graph_fd);
  //printf("closing graph_fd %x done !!\n", graph_fd);

    // lines specifying input data and output graph are declared and filled
    char line[MAX_FULL_PATH_LENGTH];

    // graph is parametrized based on qc_graph options and plotted
    //FILE* graph_fd_aux = popen("gnuplot -persist", "w");
    //FILE* graph_fd = (FILE*)calloc(1, sizeof(FILE));
    //memcpy(graph_fd, graph_fd_aux, sizeof(FILE));
    
    FILE* graph_fd = popen("gnuplot -persist", "w");
    
    if (graph_fd == NULL) {
        LOG_FATAL("Opening of file descriptor for gnuplot execution failed\n");
        return;
    }
    
    //FILE* graph_other = fopen("/tmp/prueba.txt", "w");
    //fclose(graph_other);
    //return;
printf("line: %x\n", line);
printf("graph_fd: %x\n", graph_fd);
printf("graph_fd: %i\n", graph_fd);
    sprintf(line, "set output '%s'\n", graph_filename);
printf("1 - --------------->\n");
printf("line: %s\n", line);
printf("line: %x\n", line);
printf("graph_fd: %x\n", graph_fd);
printf("graph_fd: %i\n", graph_fd);
printf("graph_fd is NULL: %i\n", (graph_fd == NULL) ? 1:0);
    fprintf(graph_fd, " ");
printf("2 - --------------->\n");            
    fprintf(graph_fd, line);
printf("3 - --------------->\n");
    fprintf(graph_fd, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
    sprintf(line, "set ylabel '%s'\n", qc_graph.ylabel);
    fprintf(graph_fd, line);
    sprintf(line, "set xlabel '%s'\n", qc_graph.xlabel);
    fprintf(graph_fd, line);
    fprintf(graph_fd, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
    sprintf(line, "set title '%s'\n", qc_graph.title);
    fprintf(graph_fd, line);

    if (qc_graph.x_autoscale == 1) {
        fprintf(graph_fd, "set autoscale x\n");
    } else {
        sprintf(line, "set xrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
        fprintf(graph_fd, line);
    }

    if (qc_graph.y_autoscale == 1) {
        fprintf(graph_fd, "set autoscale y\n");
    } else {
        sprintf(line, "set yrange [ %i : %i ] noreverse nowriteback\n", qc_graph.x_start, qc_graph.x_end);
        fprintf(graph_fd, line);
    }

    sprintf(line, "set lmargin '%i'\n", qc_graph.lmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set rmargin '%i'\n", qc_graph.rmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set tmargin '%i'\n", qc_graph.tmargin);
    fprintf(graph_fd, line);
    sprintf(line, "set bmargin '%i'\n", qc_graph.bmargin);
    fprintf(graph_fd, line);

    sprintf(line, "plot ");

    for (int i = 0; i < qc_graph.num_y_columns; i++) {
        sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", line, (i == 0 ? "" : ", "), data_filename, qc_graph.x_column, qc_graph.y_columns[i], qc_graph.y_titles[i], qc_graph.type);
    }
    fprintf(graph_fd, line);

    pclose(graph_fd);    
}

void generate_html_file_with_images(qc_report_t qc_report, int kmers_on, char* html_filename, char* report_directory, char* in_shortname) {
    char line[MAX_FULL_PATH_LENGTH];

    FILE* fd = fopen(html_filename, "w");

    fprintf(fd, "<HTML>\n");
    fprintf(fd, "<HEAD>\n");
    fprintf(fd, "<TITLE>Quality Control Results</TITLE>\n");
    fprintf(fd, "</HEAD>\n");
    fprintf(fd, "<BODY>\n");
    fprintf(fd, "<TABLE WIDTH='70%'>\n");
    fprintf(fd, "<TR><TD ALIGN='center'><U><B>QUALITY CONTROL RESULTS</B></U></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");

    // print general statistics from file
    struct stat buf;

    sprintf(line, "%s/%s.qc", report_directory, in_shortname);
    stat(line, &buf);
    char* buf_p = (char *) malloc(buf.st_size * sizeof(char));
    FILE* aux_fd = fopen(line, "r");
    fread(buf_p, 1, buf.st_size, aux_fd);
    fclose(aux_fd);

    fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");
    fprintf(fd, buf_p);
    fprintf(fd, "</PRE></TD></TR>\n");

    // print graphics
    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.nt_quality.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.quality_read.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.base_sequence.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.base_gc.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.gc_histogram.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.n_base.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.read_length.png'</TD></TR>\n", report_directory, in_shortname);
    fprintf(fd, line);

    if (kmers_on) {
        sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.kmers_positions.png'</TD></TR>\n", report_directory, in_shortname);
        fprintf(fd, line);

        // print firsts kmers counts
        sprintf(line, "%s/%s.kmers.dat", report_directory, in_shortname);

        fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");

        char buffer[1024];
        aux_fd = fopen(line, "r");
        int i = 0;
        while (i < 20 && fgets(buffer, 1024, aux_fd) != NULL) {
            fprintf(fd, buffer);
            i++;
        }
        fclose(aux_fd);

        fprintf(fd, "</PRE></TD></TR>\n");
    }


    fprintf(fd, "</TABLE>\n");
    fprintf(fd, "</BODY>\n");
    fprintf(fd, "</HTML>\n");

    fclose(fd);
}

void generate_html_file_valid_with_images(qc_report_t qc_report, int kmers_on, int cg_flag, char* html_filename, char* report_directory, char* in_shortname, int valid) {
    char line[MAX_FULL_PATH_LENGTH];
    char* str_valid_suffix = "";

    str_valid_suffix = (valid) ? (char*) ".valid" : (char*) ".invalid";

    FILE* fd = fopen(html_filename, "w");

    fprintf(fd, "<HTML>\n");
    fprintf(fd, "<HEAD>\n");
    fprintf(fd, "<TITLE>Quality Control Results</TITLE>\n");
    fprintf(fd, "</HEAD>\n");
    fprintf(fd, "<BODY>\n");
    fprintf(fd, "<TABLE WIDTH='70%'>\n");
    fprintf(fd, "<TR><TD ALIGN='center'><U><B>QUALITY CONTROL RESULTS</B></U></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");
    fprintf(fd, "<TR><TD ALIGN='center'></TD></TR>\n");

    // print general statistics from file
    struct stat buf;

    sprintf(line, "%s/%s%s%s", report_directory, in_shortname, QC_SUFFIX, str_valid_suffix);
    stat(line, &buf);

    char* buf_p = (char *) calloc(1, buf.st_size * sizeof(char));

    FILE* aux_fd = fopen(line, "r");

    fread(buf_p, 1, buf.st_size, aux_fd);
    fclose(aux_fd);

    fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");
    fprintf(fd, buf_p);
    fprintf(fd, "</PRE></TD></TR>\n");
    fprintf(fd, "<TR><TD></TD></TR>\n");
    free(buf_p);

    // print chaos game data from file
    if (cg_flag) {
        sprintf(line, "%s/%s%s%s", report_directory, in_shortname, CG_SUFFIX, str_valid_suffix);
        stat(line, &buf);

        buf_p = (char*) calloc(1, buf.st_size * sizeof(char));

        aux_fd = fopen(line, "r");

        fread(buf_p, 1, buf.st_size, aux_fd);
        fclose(aux_fd);

        fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");
        fprintf(fd, buf_p);
        fprintf(fd, "</PRE></TD></TR>\n");
        fprintf(fd, "<TR><TD></TD></TR>\n");  //Rows of separation
        fprintf(fd, "<TR><TD></TD></TR>\n");
        free(buf_p);

        // printg chaos game images
        sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s'</TD></TR>\n", qc_report.pgm_fastq_filename);
        fprintf(fd, line);

        sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s'</TD></TR>\n", qc_report.pgm_quality_filename);
        fprintf(fd, line);

        sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s'</TD></TR>\n", qc_report.pgm_diff_filename);
        fprintf(fd, line);
    }

    // print graphics
    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.nt_quality%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.quality_read%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.base_sequence%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.base_gc%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.gc_histogram%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.n_base%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.read_length%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
    fprintf(fd, line);

    if (kmers_on) {
        sprintf(line, "<TR><TD ALIGN='center'><IMG src='%s/%s.kmers_positions%s.png'</TD></TR>\n", report_directory, in_shortname, str_valid_suffix);
        fprintf(fd, line);

        // print firsts kmers counts
        sprintf(line, "%s/%s.kmers.dat%s", report_directory, in_shortname, str_valid_suffix);

        fprintf(fd, "<TR><TD ALIGN='center'><PRE>\n");

        char buffer[1024];
        aux_fd = fopen(line, "r");

        int i = 0;
        while (i < 20 && fgets(buffer, 1024, aux_fd) != NULL) {
            fprintf(fd, buffer);
            i++;
        }
        fclose(aux_fd);
        fprintf(fd, "</PRE></TD></TR>\n");
    }

    fprintf(fd, "</TABLE>\n");
    fprintf(fd, "</BODY>\n");
    fprintf(fd, "</HTML>\n");

    fclose(fd);
}
