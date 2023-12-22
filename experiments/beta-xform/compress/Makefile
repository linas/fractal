
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: remap skew \
	sampler \
	nilpot \
	moments \
	zero tribo midmap \
	fpcoef poincare \
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
nilpot.o: nilpot.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

poincare.o: poincare.C

remap: remap.o
skew: skew.o
sampler: sampler.o
nilpot: nilpot.o
moments: moments.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o

poincare: poincare.o $(BRAT)
momenta: momenta.o $(BRAT)
