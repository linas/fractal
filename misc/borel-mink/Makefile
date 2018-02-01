
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: borel

borel: borel.o stern.o


.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o: stern.h
	cc -c $^  -Wall -O2 
