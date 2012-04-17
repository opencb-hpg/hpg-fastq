

BIN = ./bin

COMMONS_LIB = ../../commons
COMMONS_CUDA_LIB = ../../commons-cuda
CONTAINERS_LIB=../../containers
FASTQ_LIB = ../../bioinfo-data/fastq

ALL = hpg-fastq

CC = gcc
CFLAGS = -Wall -O3 -std=c99 -fopenmp
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

CINCLUDES = -I. -I/opt/cuda/include -I../../bioinfo-data/fastq -I../../commons -I../../commons-cuda -I../../containers
CUINCLUDES = -I. -I../../bioinfo-data/fastq -I../../commons -I../../commons-cuda -I../../containers

NVCC = nvcc
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12

all: $(ALL)

hpg-fastq: fastq_hpc_main.o prepro.o prepro_kernel_cuda.o cuda_commons.o hpg-fastq-objects
	$(NVCC) $(NVCCFLAGS) cuda_commons.o fastq_hpc_main.o string_utils.o file_utils.o system_utils.o log.o list.o prepro_kernel_cuda.o prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o \
		fastq_batch.o fastq_batch_list.o fastq_batch_reader.o qc_batch.o prepro_batch.o result_batch.o prepro.o chaos_game.o -o $(BIN)/hpg-fastq

file_utils.o: $(COMMONS_LIB)/file_utils.h
	$(CC) $(CFLAGS) -c $(COMMONS_LIB)/file_utils.c

system_utils.o: $(COMMONS_LIB)/system_utils.h
	$(CC) $(CFLAGS) -c $(COMMONS_LIB)/system_utils.c

string_utils.o: $(COMMONS_LIB)/string_utils.h $(COMMONS_LIB)/string_utils.c
	$(CC) $(CFLAGS) -c $(COMMONS_LIB)/string_utils.c

log.o: $(COMMONS_LIB)/log.h $(COMMONS_LIB)/log.c string_utils.o
	$(CC) $(CFLAGS) -c $(COMMONS_LIB)/log.c

prepro.o: prepro.cu prepro.h *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -c prepro.cu

prepro_kernel_cuda.o: prepro_kernel_cuda.cu *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -c prepro_kernel_cuda.cu

cuda_commons.o: $(COMMONS_CUDA_LIB)/cuda_commons.cu $(COMMONS_CUDA_LIB)/cuda_commons.h *.h
	$(NVCC) $(NVCCFLAGS) -c $(COMMONS_CUDA_LIB)/cuda_commons.cu

fastq_hpc_main.o: fastq_hpc_main.c *.h
	$(CC) $(CFLAGS) $(CINCLUDES) -c fastq_hpc_main.c

hpg-fastq-objects:
	$(CC) $(CFLAGS) $(CINCLUDES) -c $(CONTAINERS_LIB)/list.c $(FASTQ_LIB)/*.c $(COMMONS_LIB)/file_utils.c $(COMMONS_LIB)/log.c  $(COMMONS_LIB)/string_utils.c \
$(COMMONS_LIB)/system_utils.c qc_batch.c qc_report.c prepro_batch.c prepro_kernel_omp.c result_batch.c chaos_game.c

clean:
	rm -f *~ \#*\# *.o $(BIN)/hpg-fastq
