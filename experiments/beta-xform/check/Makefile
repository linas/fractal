
CC = cc -Wall -g -O2 $(INCLUDES)

all: recheck force

.c.o:
	$(CC) -c $<

.o:
	$(CC) -o $@ $^ -lm -lstdc++


recheck: recheck.o
force: force.o
