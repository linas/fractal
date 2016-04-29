
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

all: clifford

clifford: clifford.o
	cc -o $@ $^ -lm

#	g++ -o $@ $< -L$(LIB)  -lfunc -lm

.c.o:
	cc -std=gnu99 -c $^ -I$(INC) 
