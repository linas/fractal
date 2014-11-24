
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

all: bessy eigenwave

bessy: schroed.o bessy.o
eigenwave: schroed.o eigenwave.o

.o:
	cc -o $* $^ -L$(LIB) -lfunc -lm -lstdc++

.cc.o:
	cc -c -O2 -I $(INC) $<
