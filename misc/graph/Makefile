
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

all: bessy eigenwave eigenzee potential

bessy: schroed.o bessy.o
eigenwave: schroed.o eigenwave.o
eigenzee: schroed.o eigenzee.o
potential: potential.o

.o:
	cc -o $* $^ -L$(LIB) -lfunc -lm -lstdc++

.cc.o:
	cc -c -O2 -I $(INC) $<
