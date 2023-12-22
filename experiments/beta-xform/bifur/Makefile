
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: bernie bigbern remap unbernie skew fp-diagram \
	carry-sum rexfer sampler \
	nilpot bitshow complex moments zero tribo midmap \
	poincare momenta

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

bernie.o: bernie.C $(FUNCDIR)/bitops.h
bigbern.o: bigbern.C
remap.o: remap.c
skew.o: skew.c
unbernie.o: unbernie.C
fp-diagram.o: fp-diagram.C
carry-sum.o: carry-sum.C
rexfer.o: rexfer.c
sampler.o: sampler.c
nilpot.o: nilpot.c
bitshow.o: bitshow.c $(FUNCDIR)/bitops.h
complex.o: complex.c
moments.o: moments.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c

poincare.o: poincare.C

bernie: bernie.o $(BRAT) $(FUNC)
bigbern: bigbern.o $(BRAT)
remap: remap.o
skew: skew.o
unbernie: unbernie.o $(BRAT)
fp-diagram: fp-diagram.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
rexfer: rexfer.o
sampler: sampler.o
nilpot: nilpot.o
bitshow: bitshow.o $(FUNC)
complex: complex.o
moments: moments.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o

poincare: poincare.o $(BRAT)
momenta: momenta.o $(BRAT)
