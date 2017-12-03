
CC = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: circle-map

.C.o:
	$(CC) -c $< -I $(GENDIR)

.o:
	$(CC) -o $@ $^ -lpthread -lm -lstdc++

circle-map.o: circle-map.C

circle-map: circle-map.o $(BRAT)
