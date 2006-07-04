############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC)


EXES= angular buddy hyperg mertens ndist series von_mangoldt

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf num.tgz *.c *.h *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $^ -lm


angular:	angular.o $(LIB)/cache.o $(LIB)/gcf.o $(LIB)/modular.o $(LIB)/moebius.o $(LIB)/totient.o
buddy:	buddy.o
hyperg:	hyperg.o $(LIB)/modular.o
mertens:	mertens.o $(LIB)/cache.o $(LIB)/gcf.o $(LIB)/moebius.o
ndist:	ndist.o $(LIB)/gcf.o $(LIB)/totient.o
series:	series.o $(LIB)/cache.o $(LIB)/gcf.o $(LIB)/modular.o $(LIB)/moebius.o $(LIB)/totient.o 
von_mangoldt:	von_mangoldt.o $(LIB)/cache.o $(LIB)/gcf.o $(LIB)/moebius.o

ndist.o:	ndist.c $(INC)/gcf.h $(INC)/totient.h 

