
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: green sho

green:	green.o polynomial.o

sho:	sho.o

green.o:	green.c
polynomial.o:	polynomial.c gaussian.h
sho.o: sho.c

.c.o:
	cc -c -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -L$(LIB) -lfunc -lgsl -lcblas -lm

