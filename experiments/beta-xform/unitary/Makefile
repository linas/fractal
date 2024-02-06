
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

PERIODIC = ../periodic

all: bigun bigfourier poly-four complex-root

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(PERIODIC)

.c.o:
	$(CC) -c $< -I $(PERIODIC)

.o:
	$(CC) -o $@ $^ -lanant -ldb -lgmp -lm -lstdc++

bigun.o: bigun.c bigseries.c bigzero.c
bigfourier.o: bigfourier.C bigseries.c
poly-four.o: poly-four.C
complex-root.o: complex-root.c bigseries.c bigzero.c

bigfourier: bigfourier.o $(BRAT)
poly-four: poly-four.o $(BRAT)
