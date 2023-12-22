
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: remap reigen skew fp-diagram \
	carry-sum rexfer treefn treestep sampler julia julie \
	nilpot bitshow complex \
	moments \
	zero tribo midmap \
	fpcoef poincare beta witt \
	momenta 

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

remap.o: remap.c
reigen.o: reigen.c
skew.o: skew.c
fp-diagram.o: fp-diagram.C
carry-sum.o: carry-sum.C
rexfer.o: rexfer.c
treefn.o: treefn.c
treestep.o: treestep.C
sampler.o: sampler.c
julia.o: julia.c
julie.o: julie.C
nilpot.o: nilpot.c
bitshow.o: bitshow.c $(FUNCDIR)/bitops.h
complex.o: complex.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

poincare.o: poincare.C
beta.o: beta.c
witt.o: witt.c

remap: remap.o
reigen: reigen.o
skew: skew.o
fp-diagram: fp-diagram.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
rexfer: rexfer.o
treefn: treefn.o
treestep: treestep.o  $(BRAT)
sampler: sampler.o
julia: julia.o
julie: julie.o $(BRAT)
nilpot: nilpot.o
bitshow: bitshow.o $(FUNC)
complex: complex.o
moments: moments.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o
fpcoef: fpcoef.o
	$(CC) -o $@ $^ -lgmp -lgomp -lm

poincare: poincare.o $(BRAT)
beta: beta.o
witt: witt.o
momenta: momenta.o $(BRAT)
