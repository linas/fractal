
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../../tools/inc
FUNC=../../../tools/lib/libfunc.a

all: reigen fp-diagram rexfer complex moments beta witt

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

reigen.o: reigen.c
fp-diagram.o: fp-diagram.C
rexfer.o: rexfer.c
complex.o: complex.c
moments.o: moments.c
beta.o: beta.c
witt.o: witt.c

reigen: reigen.o
fp-diagram: fp-diagram.o $(BRAT)
rexfer: rexfer.o
complex: complex.o
moments: moments.o
beta: beta.o
witt: witt.o
