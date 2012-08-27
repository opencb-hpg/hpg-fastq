
BIN = "./bin"

COMMONS_LIB = ./libs/commons
COMMONS_CUDA_LIB = ./libs/commons-cuda
CONTAINERS_LIB = ./libs/containers
FASTQ_LIB = ./libs/bioformats/fastq

NVCC_DISCOVER = $(shell expr `which nvcc | wc -l` \> 0)
#CUDA_LIB_PATH = $(shell echo -L$LD_LIBRARY_PATH"lib" | sed 's/:/ -L/g')

ALL = hpg-fastq

CC = gcc
CFLAGS = -Wall -O3 -std=c99 -fopenmp
#CFLAGS = -Wall -O0 -g -std=c99 -fopenmp
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

CINCLUDES = -I. -I/opt/cuda/include -I$(FASTQ_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)
CUINCLUDES = -I. -I$(FASTQ_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)

NVCC = nvcc

NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
#NVCCFLAGS = -O -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12

all: $(ALL)


ifeq ($(NVCC_DISCOVER), 1)
CUDA_OBJECTS = prepro_kernel_cuda.o cuda_commons.o
endif

ifeq ($(NVCC_DISCOVER), 1)
hpg-fastq: hpg-fastq-objects system_utils.o prepro.o fastq_hpc_main.o $(CUDA_OBJECTS)
	$(NVCC) $(NVCCFLAGS) fastq_hpc_main.o string_utils.o file_utils.o system_utils.o log.o list.o prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o \
		fastq_batch.o fastq_batch_list.o fastq_batch_reader.o qc_batch.o prepro_batch.o prepro_commons.o prepro.o chaos_game.o $(CUDA_OBJECTS) -o $(BIN)/hpg-fastq
else
hpg-fastq: hpg-fastq-objects system_utils.o prepro.o fastq_hpc_main.o
	$(CC) $(CFLAGS) fastq_hpc_main.o string_utils.o file_utils.o system_utils.o log.o list.o prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o \
		fastq_batch.o fastq_batch_list.o fastq_batch_reader.o qc_batch.o prepro_batch.o prepro_commons.o prepro.o chaos_game.o -o $(BIN)/hpg-fastq -lm
endif

ifeq ($(NVCC_DISCOVER), 1)
prepro.o: prepro.cu prepro.h *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -DCUDA_VERSION -c prepro.cu

prepro_kernel_cuda.o: prepro_kernel_cuda.cu *.h
	$(NVCC) $(NVCCFLAGS) $(CUINCLUDES) -DCUDA_VERSION -c prepro_kernel_cuda.cu

cuda_commons.o: #$(COMMONS_CUDA_LIB)/cuda_commons.cu $(COMMONS_CUDA_LIB)/cuda_commons.h *.h
	$(NVCC) $(NVCCFLAGS) -DCUDA_VERSION -c $(COMMONS_CUDA_LIB)/cuda_commons.cu
else
prepro.o: prepro.c prepro.h *.h
	$(CC) $(CFLAGS) $(CINCLUDES) -c prepro.c
endif

ifeq ($(NVCC_DISCOVER), 1)
system_utils.o: $(COMMONS_LIB)/system_utils.h
	@echo $(NVCC_DISCOVER)
	$(CC) $(CFLAGS) $(CINCLUDES) -DCUDA_VERSION -c $(COMMONS_LIB)/system_utils.c
else
system_utils.o: $(COMMONS_LIB)/system_utils.h
	@echo $(NVCC_DISCOVER)
	$(CC) $(CFLAGS) $(CINCLUDES) -c $(COMMONS_LIB)/system_utils.c
endif


main_hpg_fastq.o: main_hpg_fastq.c *.h
	$(CC) $(CFLAGS) $(CINCLUDES) -c main_hpg_fastq.c

hpg-fastq-objects:
	$(CC) $(CFLAGS) $(CINCLUDES) -c $(CONTAINERS_LIB)/list.c $(FASTQ_LIB)/*.c $(COMMONS_LIB)/file_utils.c $(COMMONS_LIB)/log.c  $(COMMONS_LIB)/string_utils.c \
qc_batch.c qc_report.c prepro_batch.c prepro_kernel_omp.c prepro_commons.c chaos_game.c

clean:
	rm -f *~ \#*\# *.o $(BIN)/hpg-fastq
