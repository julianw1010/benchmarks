###############################################################################
# Makefile to build the HashJoin workload
#
# Paper: Mitosis - Mitosis: Transparently Self-Replicating Page-Tables
#                  for Large-Memory Machines
# Authors: Reto Achermann, Jayneel Gandhi, Timothy Roscoe,
#          Abhishek Bhattacharjee, and Ashish Panwar
###############################################################################

CC = gcc

# cflags used
CFLAGS=-O3 -march=native

# common libraries used
COMMON_LIBS=-lrt -ldl -lpthread -lm

# sources
SRC=src/hashjoin.c src/murmurhash.c
DEPS=src/config.h src/murmur3.h

# the targets
all : hashjoin

###############################################################################
# HashJoin
###############################################################################

hashjoin: $(SRC) $(DEPS)
	$(CC) $(INCLUDES) $(CFLAGS) -fopenmp $(SRC) -o $@ $(COMMON_LIBS)

###############################################################################
# Clean
###############################################################################

clean:
	rm -f hashjoin
