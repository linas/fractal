
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../../tools/inc
FUNC=../../../tools/lib/libfunc.a

all: nilpot irred zero tribo midmap irred-viz irred-fraction \
	irred-gf irred-disk irred-midpoint irred-allowed irred-tacky \
	irred-height irred-tree irred-seek eventually bracket-finite \
	comb-limit slice bracket-steps good-rational

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

nilpot.o: nilpot.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

irred.o: irred.c irred-gold.c
irred-viz.o: irred-viz.c irred-gold.c
irred-fraction.o: irred-fraction.c irred-gold.c selfie.c
irred-gf.o: irred-gf.C irred-gold.c
irred-disk.o: irred-disk.C selfie.c
irred-midpoint.o: irred-midpoint.c irred-gold.c
irred-allowed.o: irred-allowed.c irred-gold.c
irred-tacky.o: irred-tacky.c irred-gold.c
irred-height.o: irred-height.c irred-gold.c
irred-tree.o: irred-tree.c selfie.c selfie-util.c selfie-tree.c
irred-seek.o: irred-seek.c irred-gold.c
eventually.o: eventually.c selfie.c selfie-rational.c selfie-tree.c selfie-util.c
bracket-finite.o: bracket-finite.c selfie.c selfie-rational.c selfie-tree.c selfie-util.c
comb-limit.o: comb-limit.c
slice.o: slice.c
bracket-steps.o: bracket-steps.c selfie.c
good-rational.o: good-rational.c selfie.c selfie-rational.c selfie-tree.c selfie-util.c

zero: zero.o $(BRAT)

irred-gf: irred-gf.o $(BRAT)
irred-disk: irred-disk.o $(BRAT)
comb-limit: comb-limit.o $(FUNC)
