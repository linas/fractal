
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../../tools/inc
FUNC=../../../tools/lib/libfunc.a

all: treefn treestep julia julie

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

treefn.o: treefn.c
treestep.o: treestep.C
julia.o: julia.c
julie.o: julie.C

treefn: treefn.o
treestep: treestep.o  $(BRAT)
julia: julia.o
julie: julie.o $(BRAT)
