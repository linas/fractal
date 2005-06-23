

#
# Linas Vepstas January 1994

TOP = /usr/home/src/tools

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/fcc -f2 -O -I $(INC) -DAIX221
CC = /usr/bin/cc -O2 -I $(INC) -DANSI_C

all: inc fileio.o

inc:	fileio.h
	touch inc
	cp fileio.h $(INC)/fileio.h

fileio.o:	fileio.c fileio.h
	$(CC) -c fileio.c
	cp fileio.o $(LIB)/fileio.o
	rm fileio.o
	touch fileio.o


clean:
	rm -f *.o junk tmp core a.out

realclean: clean
	rm -f inc


