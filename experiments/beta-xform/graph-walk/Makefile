
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

PERIODIC = ../periodic

all: walker

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(PERIODIC)

.c.o:
	$(CC) -c $< -I $(PERIODIC)

.o:
	$(CC) -o $@ $^ -lanant -ldb -lgmp -lm -lstdc++

walker.o: walker.c
