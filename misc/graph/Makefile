
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

all: schroed

schroed: schroed.o


.o:
	cc -o $* $< -L$(LIB) -lfunc -lm -lstdc++

.cc.o:
	cc -c -O2 -I $(INC) $<
