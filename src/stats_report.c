/*
 * stats_report.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "stats_report.h"


//--------------------------------------------------------------------

void report_summary(char *in_filename, stats_counters_t *counters, char *out_dir);
void report_length(char *in_filename, stats_counters_t *counters, char *out_dir);
void report_nt_content(char *in_filename, stats_counters_t *counters, char *out_dir);
void report_quality(char *in_filename, stats_counters_t *counters, char *out_dir);


//--------------------------------------------------------------------
// report graph
//--------------------------------------------------------------------

void init_report_graph(report_graph_t* graph);
void generate_gnuplot_image(report_graph_t *graph, char *data_filename, char *prefix);

//--------------------------------------------------------------------
// stats report
//--------------------------------------------------------------------

void stats_report(char *in_filename, stats_counters_t *counters, char *out_dir) {
  report_summary(in_filename, counters, out_dir);
  report_length(in_filename, counters, out_dir);
  report_quality(in_filename, counters, out_dir);
  report_nt_content(in_filename, counters, out_dir);
}

//--------------------------------------------------------------------

void report_summary(char *in_filename, stats_counters_t *counters, char *out_dir) {
  
  char path[strlen(in_filename) + strlen(out_dir) + 100];
  sprintf(path, "%s/%s.summary.txt", out_dir, in_filename);

  FILE *f = fopen(path, "w");

  fprintf(f, "-----------------------------------\n");
  fprintf(f, "      FastQ quality report\n");
  fprintf(f, "-----------------------------------\n");

  size_t num_nucleotides = counters->num_As + counters->num_Cs + 
    counters->num_Gs + counters->num_Ts + counters->num_Ns;

  fprintf(f, "FastQ filename: %s\n", in_filename);
  fprintf(f, "\n");
  fprintf(f, "Number of reads: %lu\n", counters->num_reads);
  fprintf(f, "Read length (min., mean, max.): (%i, %0.2f, %i)\n",
	 counters->min_length, 1.0f * counters->acc_length / counters->num_reads,
	 counters->max_length);
  fprintf(f, "\n");
  fprintf(f, "Mean quality = %0.2f\n", 1.0f * counters->acc_quality / counters->num_reads);
  fprintf(f, "\n");
  fprintf(f, "Nucleotide content (A, C, G, T, N)\n");
  fprintf(f, "\tA: %0.2f %\n", 100.0f * counters->num_As / num_nucleotides);
  fprintf(f, "\tT: %0.2f %\n", 100.0f * counters->num_Ts / num_nucleotides);
  fprintf(f, "\tG: %0.2f %\n", 100.0f * counters->num_Gs / num_nucleotides);
  fprintf(f, "\tC: %0.2f %\n", 100.0f * counters->num_Cs / num_nucleotides);
  fprintf(f, "\tN: %0.2f %\n", 100.0f * counters->num_Ns / num_nucleotides);
  fprintf(f, "GC content\n");
  fprintf(f, "\tCG: %0.2f %\n", 100.0f * (counters->num_Gs + counters->num_Cs) / num_nucleotides);
  fprintf(f, "\n");

  fprintf(f, "Mean quality per nucleotide position\n");
  khiter_t k;
  khash_t(32) *hacc_quals, *hcount_quals;
  size_t acc_quals, count_quals;
  hacc_quals = counters->kh_acc_quality_per_nt;
  hcount_quals = counters->kh_count_quality_per_nt;

  for (k = 0; k < counters->max_length; k++) {
    acc_quals = (kh_size(hacc_quals) && kh_exist(hacc_quals, k) ? kh_value(hacc_quals, k) : 0);
    count_quals = (kh_size(hcount_quals) && kh_exist(hcount_quals, k) ? kh_value(hcount_quals, k) : 0);

    fprintf(f, "\tpos. %i: %i\t", k + 1, acc_quals / count_quals);
    if ((k+1) % 5 == 0) fprintf(f, "\n");
  }
  fprintf(f, "\n");
  
  fclose(f);
}

//--------------------------------------------------------------------

void report_length(char *in_filename, stats_counters_t *counters, char *out_dir) {
  int name_length = strlen(in_filename) + strlen(out_dir) + 100;

  char data_filename[name_length];
  char img_prefix[name_length];

  // data
  sprintf(data_filename, "%s/%s.length.histogram.data", out_dir, in_filename);

  FILE *f = fopen(data_filename, "w");

  size_t lengths[counters->max_length + 1];
  memset(lengths, 0, (counters->max_length + 1) * sizeof(size_t));

  khiter_t k;
  khash_t(32) *h = counters->kh_length_histogram;
  for (k = kh_begin(h); k != kh_end(h); ++k) { 
    if (kh_exist(h, k)) {
      lengths[kh_key(h, k)] = kh_value(h, k);
    }
  }
  for (int i = 1 ; i <= counters->max_length; i++) {
    fprintf(f, "%i\t%i\n", i, lengths[i]);
  }

  fclose(f);

  // image
  sprintf(img_prefix, "%s/%s.length.histogram", out_dir, in_filename);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Read Length Histogram");
  strcpy(graph.xlabel, "Read length");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = counters->max_length + 1;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_nt_content(char *in_filename, stats_counters_t *counters, char *out_dir) {
  int name_length = strlen(in_filename) + strlen(out_dir) + 100;

  char data_filename[name_length];
  char img_prefix[name_length];

  // GC histogram data
  sprintf(data_filename, "%s/%s.GC.histogram.data", out_dir, in_filename);

  FILE *f = fopen(data_filename, "w");

  size_t num_cols = 100;
  size_t gc[num_cols];
  memset(gc, 0, num_cols * sizeof(size_t));

  khiter_t k;
  khash_t(32) *h = counters->kh_gc_histogram;
  for (k = kh_begin(h); k != kh_end(h); ++k) { 
    if (kh_exist(h, k)) {
      gc[kh_key(h, k)] = kh_value(h, k);
    }
  }
  for (int i = 1 ; i < num_cols; i++) {
    if (gc[i]) {
      fprintf(f, "%i\t%i\n", i, gc[i]);
    }
  }

  fclose(f);

  // GC histogram image
  sprintf(img_prefix, "%s/%s.GC.histogram", out_dir, in_filename);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "GC Content Histogram");
  strcpy(graph.xlabel, "GC content (%)");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // nucleotide content A, T, G, C, N
  khash_t(32) *hA, *hT, *hG, *hC, *hN;
  size_t num_total, num_As, num_Ts, num_Gs, num_Cs, num_Ns;
  hA = counters->kh_num_As_per_nt;
  hT = counters->kh_num_Ts_per_nt;
  hG = counters->kh_num_Gs_per_nt;
  hC = counters->kh_num_Cs_per_nt;
  hN = counters->kh_num_Ns_per_nt;

  // GC content per nucleotide position data
  sprintf(data_filename, "%s/%s.GC.per.nt.data", out_dir, in_filename);

  f = fopen(data_filename, "w");
  float float_value;
  for (k = 0; k < counters->max_length; k++) {
    num_As = (kh_size(hA) && kh_exist(hA, k) ? kh_value(hA, k) : 0);
    num_Ts = (kh_size(hT) && kh_exist(hT, k) ? kh_value(hT, k) : 0);
    num_Gs = (kh_size(hG) && kh_exist(hG, k) ? kh_value(hG, k) : 0);
    num_Cs = (kh_size(hC) && kh_exist(hC, k) ? kh_value(hC, k) : 0);
    num_Ns = (kh_size(hN) && kh_exist(hN, k) ? kh_value(hN, k) : 0);
    num_total = num_As + num_Ts + num_Gs + num_Cs + num_Ns;
    
    if ((float_value = 100.0f * (num_Gs + num_Cs) / num_total) > 1.0f) {
      fprintf(f, "%i\t%0.2f\n", k + 1, float_value);
    }
    //    printf("\tposition %i: A: %0.2f %, T: %0.2f %, G: %0.2f %, C: %0.2f %, N: %0.2f, quality = %0.2f\n", 
    //	   k, 100.0f * num_As / num_total, 100.0f * num_Ts / num_total, 100.0f * num_Gs / num_total, 
    //	   100.0f * num_Cs / num_total, 100.0f * num_Ns / num_total, 1.0f * acc_quals / count_quals);
  }

  fclose(f);

  // GC per nt position image
  sprintf(img_prefix, "%s/%s.GC.per.nt", out_dir, in_filename);
  init_report_graph(&graph);
  
  strcpy(graph.title, "GC Content per Nucleotide Position");
  strcpy(graph.xlabel, "Nucleotide position");
  strcpy(graph.ylabel, "GC content (%)");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = counters->max_length + 1;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // GC per nt position data
  sprintf(data_filename, "%s/%s.quality.per.nt.data", out_dir, in_filename);

  f = fopen(data_filename, "w");

  khash_t(32) *hacc_quals, *hcount_quals;
  size_t acc_quals, count_quals;
  hacc_quals = counters->kh_acc_quality_per_nt;
  hcount_quals = counters->kh_count_quality_per_nt;

  for (k = 0; k < counters->max_length; k++) {
    acc_quals = (kh_size(hacc_quals) && kh_exist(hacc_quals, k) ? kh_value(hacc_quals, k) : 0);
    count_quals = (kh_size(hcount_quals) && kh_exist(hcount_quals, k) ? kh_value(hcount_quals, k) : 0);

    fprintf(f, "%i\t%0.2f\n", k, 1.0f * acc_quals / count_quals);
  }

  fclose(f);

  // GC per nt position image
  sprintf(img_prefix, "%s/%s.quality.per.nt", out_dir, in_filename);
  init_report_graph(&graph);
  
  strcpy(graph.title, "Quality per Nucleotide Position");
  strcpy(graph.xlabel, "Nucleotide position");
  strcpy(graph.ylabel, "Quality");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = counters->max_length + 1;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // Nucleotide content per nucleotide position data
  sprintf(data_filename, "%s/%s.nucleotides.data", out_dir, in_filename);

  f = fopen(data_filename, "w");
  for (k = 0; k < counters->max_length; k++) {
    num_As = (kh_size(hA) && kh_exist(hA, k) ? kh_value(hA, k) : 0);
    num_Ts = (kh_size(hT) && kh_exist(hT, k) ? kh_value(hT, k) : 0);
    num_Gs = (kh_size(hG) && kh_exist(hG, k) ? kh_value(hG, k) : 0);
    num_Cs = (kh_size(hC) && kh_exist(hC, k) ? kh_value(hC, k) : 0);
    num_Ns = (kh_size(hN) && kh_exist(hN, k) ? kh_value(hN, k) : 0);
    num_total = num_As + num_Ts + num_Gs + num_Cs + num_Ns;
    
    fprintf(f, "%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", 
	    k + 1, 100.0f * num_As / num_total, 100.0f * num_Ts / num_total, 100.0f * num_Gs / num_total, 
	    100.0f * num_Cs / num_total, 100.0f * num_Ns / num_total);
  }

  fclose(f);

  // N per nt position image
  sprintf(img_prefix, "%s/%s.nucleotides", out_dir, in_filename);
  init_report_graph(&graph);
  
  strcpy(graph.title, "Nucleotide Content per Position");
  strcpy(graph.xlabel, "Nucleotide position");
  strcpy(graph.ylabel, "Nucleotide content (%)");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = counters->max_length + 1;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 5;
  graph.y_columns[0] = 2;
  graph.y_columns[1] = 3;
  graph.y_columns[2] = 4;
  graph.y_columns[3] = 5;
  graph.y_columns[4] = 6;
  strcpy(graph.y_titles[0], "A %");  
  strcpy(graph.y_titles[1], "T %");
  strcpy(graph.y_titles[2], "G %");
  strcpy(graph.y_titles[3], "C %");
  strcpy(graph.y_titles[4], "N %");

  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_quality(char *in_filename, stats_counters_t *counters, char *out_dir) {
  int name_length = strlen(in_filename) + strlen(out_dir) + 100;

  char data_filename[name_length];
  char img_prefix[name_length];

  // quality histogram data
  sprintf(data_filename, "%s/%s.quality.histogram.data", out_dir, in_filename);

  FILE *f = fopen(data_filename, "w");

  int num_cols = 100;
  size_t quals[num_cols];
  memset(quals, 0, num_cols * sizeof(size_t));

  khiter_t k;
  khash_t(32) *h = counters->kh_quality_histogram;
  for (k = kh_begin(h); k != kh_end(h); ++k) { 
    if (kh_exist(h, k)) {
      quals[kh_key(h, k)] = kh_value(h, k);
    }
  }
  for (int i = 1 ; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i, quals[i]);
  }

  fclose(f);

  // quality histogram image
  sprintf(img_prefix, "%s/%s.quality.histogram", out_dir, in_filename);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Quality Histogram");
  strcpy(graph.xlabel, "Quality");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // quality per nt position
  sprintf(data_filename, "%s/%s.quality.per.nt.data", out_dir, in_filename);

  f = fopen(data_filename, "w");

  num_cols = counters->max_length;

  //  khiter_t k;
  khash_t(32) *hacc_quals, *hcount_quals;
  size_t acc_quals, count_quals;
  hacc_quals = counters->kh_acc_quality_per_nt;
  hcount_quals = counters->kh_count_quality_per_nt;

  for (k = 0; k < counters->max_length; k++) {
    acc_quals = (kh_size(hacc_quals) && kh_exist(hacc_quals, k) ? kh_value(hacc_quals, k) : 0);
    count_quals = (kh_size(hcount_quals) && kh_exist(hcount_quals, k) ? kh_value(hcount_quals, k) : 0);

    fprintf(f, "%i\t%0.2f\n", k, 1.0f * acc_quals / count_quals);
  }

  fclose(f);

  // quality per nucleotide image
  sprintf(img_prefix, "%s/%s.quality.per.nt", out_dir, in_filename);
  init_report_graph(&graph);
  
  strcpy(graph.title, "Quality per Nucleotide Position");
  strcpy(graph.xlabel, "Nucleotide position");
  strcpy(graph.ylabel, "Quality");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------
// report graph
//--------------------------------------------------------------------

void init_report_graph(report_graph_t* graph) {
    graph->x_autoscale = 1;
    graph->x_start = 1;
    graph->x_end = 100;
    graph->y_autoscale = 1;
    graph->y_start = 0;
    graph->y_end = 100;
    graph->lmargin = 10;
    graph->rmargin = 4;
    graph->tmargin = 3;
    graph->bmargin = 4;
    *graph->title = 0;
    *graph->xlabel = 0;
    *graph->ylabel = 0;
    *graph->type = 0;
    graph->x_column = 0;
    graph->num_y_columns = 1;
    graph->y_columns[0] = 1;
}

//--------------------------------------------------------------------

void generate_gnuplot_image(report_graph_t *graph, char *data_filename, char *prefix) {
    // lines specifying input data and output graph are declared and filled
    char line[2048];
    
    char gnuplot_filename[strlen(prefix) + 100];
    sprintf(gnuplot_filename, "%s.gnuplot", prefix);

    // open the file for writing gnuplot lines
    FILE* f = fopen(gnuplot_filename, "w");
    
    if (f == NULL) {
      LOG_FATAL("Impossible save file for BAM report");
    }

    sprintf(line, "set output '%s.png'\n", prefix);
    fprintf(f, line);
    fprintf(f, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
    sprintf(line, "set ylabel '%s'\n", graph->ylabel);
    fprintf(f, line);
    sprintf(line, "set xlabel '%s'\n", graph->xlabel);
    fprintf(f, line);
    fprintf(f, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
    sprintf(line, "set title '%s'\n", graph->title);
    fprintf(f, line);

    if (graph->x_autoscale == 1) {
      fprintf(f, "set autoscale x\n");
    } else {
      sprintf(line, "set xrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
      fprintf(f, line);
    }

    if (graph->y_autoscale == 1) {
      fprintf(f, "set autoscale y\n");
    } else {
      sprintf(line, "set yrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
      fprintf(f, line);
    }

    sprintf(line, "set lmargin '%i'\n", graph->lmargin);
    fprintf(f, line);
    sprintf(line, "set rmargin '%i'\n", graph->rmargin);
    fprintf(f, line);
    sprintf(line, "set tmargin '%i'\n", graph->tmargin);
    fprintf(f, line);
    sprintf(line, "set bmargin '%i'\n", graph->bmargin);
    fprintf(f, line);

    sprintf(line, "plot ");

    for (int i = 0; i < graph->num_y_columns; i++) {
      sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", 
	      line, (i == 0 ? "" : ", "), data_filename, 
	      graph->x_column, graph->y_columns[i], 
	      graph->y_titles[i], graph->type);
    }
    fprintf(f, line);
    fprintf(f, "\n");

    fclose(f);    
    
    // build the command line by calling gnuplot followed by is instruction file
    // execute command line: gnuplot filename.gnuplot
    sprintf(line, "gnuplot %s;", gnuplot_filename);
    system(line);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
