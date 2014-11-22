

TOP = /src/tools

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/cc -O2 -O -I $(INC) -DAIXV3
CC = /bin/fcc -f2 -O -I $(INC)

all: inc gl.o

inc:	gl.h 
	touch inc
	cp gl.h $(INC)/gl.h

gl.o:	gl.c gl.h
	$(CC) -c gl.c
	cp gl.o $(LIB)/gl.o
	rm gl.o
	touch gl.o

