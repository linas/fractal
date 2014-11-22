#
#
# HISTORY:

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/fcc -f2 -O -I $(INC)
CC = /bin/cc -O -I $(INC)

all: inc stoch.o

inc:	stoch.h
	touch inc
	cp stoch.h $(INC)/stoch.h

stoch.o:	stoch.c stoch.h
	$(CC) -c stoch.c
	cp stoch.o $(LIB)/stoch.o
	rm stoch.o
	touch stoch.o

clean:
	cd test; make clean
	cd proto; make clean
	rm -f tmp junk core a.out *.o inc

realclean:	clean
	cd test; make realclean
	cd proto; make realclean
 
