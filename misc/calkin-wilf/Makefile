
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: fusc stern

fusc: fusc.o
stern: stern.o

.o:
	cc -o $@ $^ -L $(LIB) -lfunc -lm -lstdc++

.cc.o:
	cc -c $^  -Wall -O2 -I $(INC)
