
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: bernie bigbern remap unbernie reigen skew fp-diagram \
	carry-sum rexfer treefn treestep sampler julia julie psi \
	matrix psieigen

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $<

.o:
	$(CC) -o $@ $^ $(FUNC) -lgmp -lpthread -lm -lstdc++

bernie.o: bernie.C $(FUNCDIR)/bitops.h
bigbern.o: bigbern.C
remap.o: remap.c
reigen.o: reigen.c
skew.o: skew.c
unbernie.o: unbernie.C
fp-diagram.o: fp-diagram.C
carry-sum.o: carry-sum.C
rexfer.o: rexfer.c
treefn.o: treefn.c
treestep.o: treestep.C
sampler.o: sampler.c
julia.o: julia.c
julie.o: julie.C
psi.o: psi.c
matrix.o: matrix.C psi.c
psiegen.o: psiegen.c psi.c

bernie: bernie.o $(BRAT) $(FUNC)
bigbern: bigbern.o # $(BRAT)
remap: remap.o
reigen: reigen.o
skew: skew.o
unbernie: unbernie.o $(BRAT)
fp-diagram: fp-diagram.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
rexfer: rexfer.o
treefn: treefn.o
treestep: treestep.o  $(BRAT)
sampler: sampler.o
julia: julia.o
julie: julie.o $(BRAT)
psi: psi.o
matrix: matrix.o $(BRAT)
psieigen: psieigen.o -lgsl -lgslcblas
