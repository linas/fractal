
#
#
# HISTORY:
# created Linas Vepstas January 1994

TOP=../../..

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/cc -O2 -O -I $(INC) -DAIXV3
# CC = /bin/fcc -f2 -O -I $(INC)
CC = cc -O2 -I $(INC) -DANSI_C -DLINUX

RM= rm -f

all: inc Farey.o FareyTree.o Prime.o

inc:	Farey.h FareyTree.h Prime.h
	touch inc
	cp Farey.h $(INC)/Farey.h
	cp FareyTree.h $(INC)/FareyTree.h
	cp Prime.h $(INC)/Prime.h

Farey.o:	Farey.C Farey.h
	$(CC) -c Farey.C
	cp Farey.o $(LIB)/Farey.o
	rm Farey.o
	touch Farey.o

FareyTree.o:	FareyTree.C FareyTree.h
	$(CC) -c FareyTree.C
	cp FareyTree.o $(LIB)/FareyTree.o
	rm FareyTree.o
	touch FareyTree.o

Prime.o:	Prime.c Prime.h
	$(CC) -c Prime.c
	cp Prime.o $(LIB)/Prime.o
	rm Prime.o
	touch Prime.o

clean:
	cd test; make clean
	cd maps; make clean
	cd prods; make clean
	cd loops; make clean
	$(RM) core a.out *.o tmp junk glop inc

realclean: clean
	cd test; make realclean
	cd maps; make realclean
	cd prods; make realclean
	cd loops; make realclean
	$(RM) Prime inc

