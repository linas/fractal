
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: measure growth bound upper box

measure: measure.o
growth: growth.o
bound: bound.o
upper: upper.o
box: box.o


.o:
	cc -o $@ $^ -L $(LIB) -lfunc -lm -lstdc++

.cc.o:
	cc -c $^  -Wall -O2 -I $(INC)
