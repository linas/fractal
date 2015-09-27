
TOP = ../../tools
INC = $(TOP)/inc
INCLUDES = -I ../../generate

LIB = $(TOP)/lib
HIST=../../generate/brat.o ../../generate/util.o



all: lattice

lattice: lattice.o $(HIST)


.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
