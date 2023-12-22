
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../../tools/inc
FUNC=../../../tools/lib/libfunc.a

all: bernie bigbern unbernie \
	carry-sum bitshow momenta

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

bernie.o: bernie.C $(FUNCDIR)/bitops.h
bigbern.o: bigbern.C
unbernie.o: unbernie.C
carry-sum.o: carry-sum.C
bitshow.o: bitshow.c $(FUNCDIR)/bitops.h

bernie: bernie.o $(BRAT) $(FUNC)
bigbern: bigbern.o $(BRAT)
unbernie: unbernie.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
bitshow: bitshow.o $(FUNC)

momenta: momenta.o $(BRAT)
