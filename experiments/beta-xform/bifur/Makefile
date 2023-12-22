
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: bernie bigbern unbernie skew fp-diagram \
	carry-sum rexfer sampler \
	bitshow complex moments midmap \
	momenta

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

bernie.o: bernie.C $(FUNCDIR)/bitops.h
bigbern.o: bigbern.C
skew.o: skew.c
unbernie.o: unbernie.C
fp-diagram.o: fp-diagram.C
carry-sum.o: carry-sum.C
rexfer.o: rexfer.c
sampler.o: sampler.c
bitshow.o: bitshow.c $(FUNCDIR)/bitops.h
complex.o: complex.c
moments.o: moments.c
midmap.o: midmap.c

bernie: bernie.o $(BRAT) $(FUNC)
bigbern: bigbern.o $(BRAT)
skew: skew.o
unbernie: unbernie.o $(BRAT)
fp-diagram: fp-diagram.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
rexfer: rexfer.o
sampler: sampler.o
bitshow: bitshow.o $(FUNC)
complex: complex.o
moments: moments.o
midmap: midmap.o

momenta: momenta.o $(BRAT)
