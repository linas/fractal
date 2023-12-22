
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: remap skew \
	treefn treestep sampler julia julie \
	nilpot \
	moments \
	irred zero tribo midmap \
	fpcoef poincare \
	momenta irred-viz irred-fraction \
	irred-gf irred-disk irred-midpoint irred-allowed irred-tacky \
	irred-height irred-tree irred-seek eventually

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

remap.o: remap.c
skew.o: skew.c
treefn.o: treefn.c
treestep.o: treestep.C
sampler.o: sampler.c
julia.o: julia.c
julie.o: julie.C
nilpot.o: nilpot.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

poincare.o: poincare.C
irred.o: irred.c irred-gold.c
irred-viz.o: irred-viz.c irred-gold.c
irred-fraction.o: irred-fraction.c irred-gold.c selfie.c
irred-gf.o: irred-gf.C irred-gold.c
irred-disk.o: irred-disk.C irred-gold.c
irred-midpoint.o: irred-midpoint.c irred-gold.c
irred-allowed.o: irred-allowed.c irred-gold.c
irred-tacky.o: irred-tacky.c irred-gold.c
irred-height.o: irred-height.c irred-gold.c
irred-tree.o: irred-tree.c irred-gold.c selfie.c
irred-seek.o: irred-seek.c irred-gold.c
eventually.o: eventually.c irred-gold.c

remap: remap.o
skew: skew.o
treefn: treefn.o
treestep: treestep.o  $(BRAT)
sampler: sampler.o
julia: julia.o
julie: julie.o $(BRAT)
nilpot: nilpot.o
moments: moments.o
irred: irred.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o
fpcoef: fpcoef.o
	$(CC) -o $@ $^ -lgmp -lgomp -lm

poincare: poincare.o $(BRAT)
momenta: momenta.o $(BRAT)
irred-gf: irred-gf.o $(BRAT)
irred-disk: irred-disk.o $(BRAT)
