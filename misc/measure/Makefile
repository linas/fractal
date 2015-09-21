
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: measure growth

measure: measure.o stern.o
growth: growth.o stern.o


.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o: stern.h
	cc -c $^  -Wall -O2 
