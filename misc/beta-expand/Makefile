
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a
LAPACK=/usr/lib/x86_64-linux-gnu/liblapacke.a \
       /usr/lib/x86_64-linux-gnu/liblapack.a \
       /usr/lib/x86_64-linux-gnu/blas/libblas.a \
       /usr/lib/gcc/x86_64-linux-gnu/6/libgfortran.a \
       /usr/lib/x86_64-linux-gnu/libquadmath.so.0

all: sidorov sidorov-misc sidorov-density sidorov-bifur \
     sidorov-length branching-measure extended-length extended-measure \
     ext-misc ext-iterate ext-debug extended-debug

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

sidorov.o: sidorov.C
sidorov-misc.o: sidorov-misc.C sidorov-big.C
sidorov-density.o: sidorov-density.C extended-density.C sidorov-big.C
sidorov-bifur.o: sidorov-bifur.C sidorov-big.C
branching-measure.o: extended-density.C sidorov-big.C emrun.C
extended-length.o: extended-length.C extended-density.C sidorov-big.C
extended-measure.o: extended-measure.C extended-density.C sidorov-big.C
ext-misc.o: ext-misc.C emrun.C
ext-iterate.o: ext-iterate.C emrun.C
ext-debug.o: ext-debug.C emrun.C
extended-debug.o: extended-debug.C emrun.C sidorov-big.C

sidorov: sidorov.o # $(BRAT)

sidorov-misc: sidorov-misc.o
sidorov-length: sidorov-length.o
	$(CC) -o $@ $^ -lgmpxx -lgmp -lm -lstdc++

sidorov-density: sidorov-density.o
	$(CC) -o $@ $^ $(BRAT) -lgmpxx -lgmp -lm -lpthread -lstdc++

sidorov-bifur: sidorov-bifur.o
	$(CC) -o $@ $^ $(BRAT) -lgmpxx -lgmp -lm -lpthread -lstdc++
