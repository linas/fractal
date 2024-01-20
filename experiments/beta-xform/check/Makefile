
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

PERIODIC = ../periodic

all: recheck force orbit visitation borel reforce tsum unrolling

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR)

.c.o:
	$(CC) -c $< -I $(PERIODIC)

.o:
	$(CC) -o $@ $^ -lm -lstdc++

visitation.o: visitation.c visits.c

recheck: recheck.o
force: force.o
orbit: orbit.o
visitation: visitation.o
borel: borel.o $(BRAT)
reforce: reforce.o
tsum: tsum.o
unrolling: unrolling.o
