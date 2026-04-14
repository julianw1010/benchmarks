###############################################################################
# Makefile to build the HashJoin workload
#
# Paper: Mitosis - Mitosis: Transparently Self-Replicating Page-Tables
#                  for Large-Memory Machines
# Authors: Reto Achermann, Jayneel Gandhi, Timothy Roscoe,
#          Abhishek Bhattacharjee, and Ashish Panwar
###############################################################################

CC = gcc
CXX = g++

# cflags used
CFLAGS=-O3 -march=native
CXXFLAGS=-O3 -march=native

# common libraries used
COMMON_LIBS=-lrt -ldl -lpthread -lm

# sources for single- and multi-threaded builds
SRC_M=src/hashjoin.c src/murmurhash.c
SRC_S=src/hashjoin.c src/murmurhash.c
DEPS=src/config.h src/murmur3.h

# the targets
all : bench_hashjoin_st bench_hashjoin_mt

###############################################################################
# HashJoin
###############################################################################

bench_hashjoin_st: $(SRC_S) $(DEPS)
	$(CC) $(INCLUDES) $(CFLAGS) $(SRC_S) -o $@  $(COMMON_LIBS)

bench_hashjoin_mt:  $(SRC_M) $(DEPS)
	$(CC) $(INCLUDES) $(CFLAGS) -fopenmp $(SRC_M) -o $@  $(COMMON_LIBS)

###############################################################################
# Clean
###############################################################################

clean:
	rm -f bench_hashjoin_st bench_hashjoin_mt
