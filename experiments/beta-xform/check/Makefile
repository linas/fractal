
CC = cc -Wall -g -O2 $(INCLUDES)

all: recheck

.c.o:
	$(CC) -c $<

.o:
	$(CC) -o $@ $^ -lm -lstdc++


recheck: recheck.o

