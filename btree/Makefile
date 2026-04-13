CC      = gcc
CFLAGS  = -O2 -std=c99 -fopenmp -D_POSIX_C_SOURCE=200112L
LDFLAGS = -fopenmp
TARGET  = btree
SRC     = btree.c

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)
