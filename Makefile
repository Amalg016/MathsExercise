CC = gcc
CFLAGS = -Wall -O2

all:
	$(CC) $(CFLAGS) matrix.c main.c -o matrix -lm

clean:
	rm -f matrix
