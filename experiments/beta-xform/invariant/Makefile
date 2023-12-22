
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: reigen skew fp-diagram \
	rexfer sampler \
	nilpot complex moments zero tribo midmap \
	fpcoef poincare beta witt momenta 

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

reigen.o: reigen.c
skew.o: skew.c
fp-diagram.o: fp-diagram.C
rexfer.o: rexfer.c
sampler.o: sampler.c
nilpot.o: nilpot.c
complex.o: complex.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

poincare.o: poincare.C
beta.o: beta.c
witt.o: witt.c

reigen: reigen.o
skew: skew.o
fp-diagram: fp-diagram.o $(BRAT)
rexfer: rexfer.o
sampler: sampler.o
nilpot: nilpot.o
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
