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

all: inc Farey.o FareyTree.o gcf.o modular.o moebius.o Prime.o totient.o

inc:	Farey.h FareyTree.h gcf.h modular.h moebius.h Prime.h totient.h
	touch inc
	cp Farey.h $(INC)/Farey.h
	cp FareyTree.h $(INC)/FareyTree.h
	cp gcf.h $(INC)/gcf.h
	cp moebius.h $(INC)/moebius.h
	cp modular.h $(INC)/modular.h
	cp Prime.h $(INC)/Prime.h
	cp totient.h $(INC)/totient.h

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

gcf.o:	gcf.c gcf.h
	$(CC) -c gcf.c
	cp gcf.o $(LIB)/gcf.o
	rm gcf.o
	touch gcf.o

modular.o:	modular.c modular.h
	$(CC) -c modular.c
	cp modular.o $(LIB)/modular.o
	rm modular.o
	touch modular.o

moebius.o:	moebius.c moebius.h
	$(CC) -c moebius.c
	cp moebius.o $(LIB)/moebius.o
	rm moebius.o
	touch moebius.o

Prime.o:	Prime.c Prime.h
	$(CC) -c Prime.c
	cp Prime.o $(LIB)/Prime.o
	rm Prime.o
	touch Prime.o

totient.o:	totient.c totient.h
	$(CC) -c totient.c
	cp totient.o $(LIB)/totient.o
	rm totient.o
	touch totient.o

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

dist:
	tar -zcvf x *.C *.h *.c Makefile .cvsignore
