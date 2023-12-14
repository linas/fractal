
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a
LAPACK=/usr/lib/x86_64-linux-gnu/liblapacke.a \
       /usr/lib/x86_64-linux-gnu/liblapack.a \
       /usr/lib/x86_64-linux-gnu/blas/libblas.a \
       /usr/lib/gcc/x86_64-linux-gnu/12/libgfortran.a \
       /usr/lib/x86_64-linux-gnu/libquadmath.so.0

all: minside

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

minside.o: minside.C

minside: minside.o $(BRAT)
