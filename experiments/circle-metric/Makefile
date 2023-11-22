
CC = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: circle measure

.C.o:
	$(CC) -c $< -I $(GENDIR)

.o:
	$(CC) -o $@ $^ -lpthread -lm -lstdc++

circle.o: circle.C
measure.o: measure.C

circle: circle.o $(BRAT)
measure: measure.o
