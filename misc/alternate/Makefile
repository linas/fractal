
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

all: alt-binary alt-frac

alt-binary: alt-binary.o
	cc -o $@ $^ -lm

alt-frac: alt-frac.o
	g++ -o $@ $< -L$(LIB)  -lfunc -lm

.C.o:
	cc -c $^ -I$(INC) 
