
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: lattice

lattice: lattice.o


.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o:
	cc -c $^  -Wall -O2 
