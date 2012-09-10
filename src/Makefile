CC = gcc
CFLAGS = -std=c99 -O3 -march=native -Wall -fopenmp

BIN_DIR = "../bin"

LIBS_DIR = $(PWD)/../libs
COMMONS_LIB = $(LIBS_DIR)/commons
CONTAINERS_LIB = $(LIBS_DIR)/containers
BIOFORMATS_LIB = $(LIBS_DIR)/bioformats
COMMONS_CUDA_LIB = $(LIBS_DIR)/commons-cuda

NVCC_DISCOVER = $(shell expr `which nvcc | wc -l` \> 0)
#CUDA_LIB_PATH = $(shell echo -L$LD_LIBRARY_PATH"lib" | sed 's/:/ -L/g')

ALL = hpg-fastq

#CFLAGS = -Wall -O0 -g -std=c99 -fopenmp
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

INCLUDES = -I . -I$(BIOFORMATS_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)
LIBS = -L/usr/lib/x86_64-linux-gnu -lm -fopenmp -lconfig -largtable2

## To remove ##
NVCC = nvcc
NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
CUINCLUDES = -I. -I/opt/cuda/include -I$(BIOFORMATS_LIB) -I$(COMMONS_LIB) -I$(COMMONS_CUDA_LIB) -I$(CONTAINERS_LIB)
#NVCCFLAGS = -O -Xptxas -v  -arch=sm_20 -Xcompiler -fopenmp
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12

all: $(ALL)


ifeq ($(NVCC_DISCOVER), 1)
CUDA_OBJECTS = prepro_kernel_cuda.o cuda_commons.o
endif

ifeq ($(NVCC_DISCOVER), 1)
hpg-fastq: hpg-fastq-objects system_utils.o prepro.o main_hpg_fastq.o $(CUDA_OBJECTS)
	$(NVCC) $(NVCCFLAGS) main_hpg_fastq.o string_utils.o file_utils.o system_utils.o log.o list.o \
	prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o fastq_batch.o fastq_batch_list.o \
	fastq_batch_reader.o qc_batch.o prepro_batch.o prepro_commons.o prepro.o chaos_game.o $(CUDA_OBJECTS) -o $(BIN_DIR)/hpg-fastq
else
hpg-fastq: hpg-fastq-objects system_utils.o prepro.o main_hpg_fastq.o
	$(CC) $(CFLAGS) main_hpg_fastq.o string_utils.o file_utils.o system_utils.o log.o list.o \
	prepro_kernel_omp.o qc_report.o fastq_file.o fastq_read.o fastq_batch.o fastq_batch_list.o \
	fastq_batch_reader.o qc_batch.o prepro_batch.o prepro_commons.o prepro.o chaos_game.o options.o -o $(BIN_DIR)/hpg-fastq $(LIBS)
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
	$(CC) $(CFLAGS) $(INCLUDES) -c prepro.c
endif

ifeq ($(NVCC_DISCOVER), 1)
system_utils.o: $(COMMONS_LIB)/system_utils.h
	@echo $(NVCC_DISCOVER)
	$(CC) $(CFLAGS) $(INCLUDES) -DCUDA_VERSION -c $(COMMONS_LIB)/system_utils.c
else
system_utils.o: $(COMMONS_LIB)/system_utils.*
	@echo $(NVCC_DISCOVER)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(COMMONS_LIB)/system_utils.c
endif


main_hpg_fastq.o: main_hpg_fastq.c *.h
	$(CC) $(CFLAGS) $(INCLUDES) -c main_hpg_fastq.c

hpg-fastq-objects:
	$(CC) $(CFLAGS) $(INCLUDES) -c $(CONTAINERS_LIB)/list.c $(BIOFORMATS_LIB)/fastq/*.c $(COMMONS_LIB)/file_utils.c $(COMMONS_LIB)/log.c $(COMMONS_LIB)/string_utils.c \
qc_batch.c qc_report.c prepro_batch.c prepro_kernel_omp.c prepro_commons.c chaos_game.c options.c

clean:
	rm -f *~ \#*\# *.o $(BIN)/hpg-fastq
