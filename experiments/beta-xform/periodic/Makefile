
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: skew \
	sampler \
	nilpot \
	moments \
	irred zero tribo midmap \
	momenta irred-viz irred-fraction \
	irred-gf irred-disk irred-midpoint irred-allowed irred-tacky \
	irred-height irred-tree irred-seek eventually

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

skew.o: skew.c
sampler.o: sampler.c
nilpot.o: nilpot.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

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

skew: skew.o
sampler: sampler.o
nilpot: nilpot.o
moments: moments.o
irred: irred.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o

momenta: momenta.o $(BRAT)
irred-gf: irred-gf.o $(BRAT)
irred-disk: irred-disk.o $(BRAT)
