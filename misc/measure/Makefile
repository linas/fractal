
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: measure growth bound upper

measure: measure.o stern.o
growth: growth.o stern.o
bound: bound.o stern.o
upper: upper.o stern.o


.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o: stern.h
	cc -c $^  -Wall -O2 
