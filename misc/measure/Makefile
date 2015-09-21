
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: stern

stern: stern.o

.o:
	cc -o $@ $^ -lm -lstdc++

.cc.o:
	cc -c $^  -Wall -O2 
