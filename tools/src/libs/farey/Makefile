
#
# Linas Vepstas January 1994

TOP = ../../..

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/cc -O2 -O -I $(INC) -DAIXV3
# CC = /bin/fcc -f2 -O -I $(INC)
CC = /usr/bin/cc  -O2 -I $(INC) -DANSI_C

all: inc Farey.o Prime.o

inc:	Farey.h Prime.h
	touch inc
	cp Farey.h $(INC)/Farey.h
	cp Prime.h $(INC)/Prime.h

Farey.o:	Farey.c Farey.h
	$(CC) -c Farey.c
	cp Farey.o $(LIB)/Farey.o
	rm Farey.o
	touch Farey.o

Prime.o:	Prime.c Prime.h
	$(CC) -c Prime.c
	cp Prime.o $(LIB)/Prime.o
	rm Prime.o
	touch Prime.o

clean:
	rm -f *.o core tmp a.out junk

realclean:	clean


