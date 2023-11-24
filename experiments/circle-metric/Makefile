
CC = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: circle measure topo cyli bifurcation

.C.o:
	$(CC) -c $< -I $(GENDIR)

.o:
	$(CC) -o $@ $^ -lpthread -lm -lstdc++

circle: circle.o $(BRAT)
measure: measure.o
topo: topo.o
cyli: cyli.o
bifurcation: bifurcation.o
