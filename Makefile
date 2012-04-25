

BIN = ./bin

COMMONS_LIB = ./libs/commons
COMMONS_CUDA_LIB = ./libs/commons-cuda
CONTAINERS_LIB = ./libs/containers
FASTQ_LIB = ./libs/bioformats/fastq

#NVCC_DISCOVER = $(shell which nvcc | wc -l)
NVCC_DISCOVER = $(shell expr `which nvcc | wc -l` \> 0)

ALL = hpg-fastq

CC = gcc
CFLAGS = -Wall -O3 -std=c99 -fopenmp
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

CINCLUDES = -I. -I/opt/cuda/include -I$(FASTQ_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)
CUINCLUDES = -I. -I$(FASTQ_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)

ifeq ($(NVCC_DISCOVER), 1)
	NVCC = nvcc
else
	NVCC = gcc
endif

NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12

all: $(ALL)

hpg-fastq: fastq_hpc_main.o prepro.o prepro_kernel_cuda.o cuda_commons.o hpg-fastq-objects
	$(NVCC) $(NVCCFLAGS) cuda_commons.o fastq_hpc_main.o string_utils.o file_utils.o system_utils.o log.o list.o prepro_kernel_cuda.o prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o \
		fastq_batch.o fastq_batch_list.o fastq_batch_reader.o qc_batch.o prepro_batch.o prepro.o chaos_game.o -o $(BIN)/hpg-fastq

file_utils.o: $(COMMONS_LIB)/file_utils.h
	$(CC) $(CFLAGS) -DCUDA_VERSION -c $(COMMONS_LIB)/file_utils.c

system_utils.o: $(COMMONS_LIB)/system_utils.h
	$(CC) $(CFLAGS) -DCUDA_VERSION -c $(COMMONS_LIB)/system_utils.c

string_utils.o: $(COMMONS_LIB)/string_utils.h $(COMMONS_LIB)/string_utils.c
	$(CC) $(CFLAGS) -DCUDA_VERSION -c $(COMMONS_LIB)/string_utils.c

log.o: $(COMMONS_LIB)/log.h $(COMMONS_LIB)/log.c string_utils.o
	$(CC) $(CFLAGS) -DCUDA_VERSION -c $(COMMONS_LIB)/log.c

prepro.o: prepro.cu prepro.h *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -DCUDA_VERSION -c prepro.cu

prepro_kernel_cuda.o: prepro_kernel_cuda.cu *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -DCUDA_VERSION -c prepro_kernel_cuda.cu

cuda_commons.o: $(COMMONS_CUDA_LIB)/cuda_commons.cu $(COMMONS_CUDA_LIB)/cuda_commons.h *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c $(COMMONS_CUDA_LIB)/cuda_commons.cu

fastq_hpc_main.o: fastq_hpc_main.c *.h
	$(CC) $(CFLAGS) $(CINCLUDES) -DCUDA_VERSION -c fastq_hpc_main.c

hpg-fastq-objects:
	$(CC) $(CFLAGS) $(CINCLUDES) -c $(CONTAINERS_LIB)/list.c $(FASTQ_LIB)/*.c $(COMMONS_LIB)/file_utils.c $(COMMONS_LIB)/log.c  $(COMMONS_LIB)/string_utils.c \
$(COMMONS_LIB)/system_utils.c qc_batch.c qc_report.c prepro_batch.c prepro_kernel_omp.c chaos_game.c

clean:
	rm -f *~ \#*\# *.o $(BIN)/hpg-fastq
