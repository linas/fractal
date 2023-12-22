
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: remap skew \
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

remap.o: remap.c
skew.o: skew.c
sampler.o: sampler.c
moments.o: moments.c
midmap.o: midmap.c

remap: remap.o
skew: skew.o
sampler: sampler.o
moments: moments.o
midmap: midmap.o

momenta: momenta.o $(BRAT)
