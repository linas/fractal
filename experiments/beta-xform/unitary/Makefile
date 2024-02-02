
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

PERIODIC = ../periodic

all: bigun

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR)

.c.o:
	$(CC) -c $< -I $(PERIODIC)

.o:
	$(CC) -o $@ $^ -lgmp -lm -lstdc++

bigun.o: bigun.c
