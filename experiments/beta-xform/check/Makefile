
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

all: recheck force orbit visitation borel

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR)

.c.o:
	$(CC) -c $<

.o:
	$(CC) -o $@ $^ -lm -lstdc++


recheck: recheck.o
force: force.o
orbit: orbit.o
visitation: visitation.o
borel: borel.o $(BRAT)
