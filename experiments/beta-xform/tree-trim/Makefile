
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: skew \
	treefn treestep sampler julia julie \
	moments \
	zero tribo midmap \
	momenta

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

skew.o: skew.c
treefn.o: treefn.c
treestep.o: treestep.C
sampler.o: sampler.c
julia.o: julia.c
julie.o: julie.C
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c


skew: skew.o
treefn: treefn.o
treestep: treestep.o  $(BRAT)
sampler: sampler.o
julia: julia.o
julie: julie.o $(BRAT)
moments: moments.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o

momenta: momenta.o $(BRAT)
