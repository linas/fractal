
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: skew \
	sampler \
	moments \
	midmap \
	momenta 

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

skew.o: skew.c
sampler.o: sampler.c
moments.o: moments.c
midmap.o: midmap.c

skew: skew.o
sampler: sampler.o
moments: moments.o

midmap: midmap.o

poincare: poincare.o $(BRAT)
momenta: momenta.o $(BRAT)
