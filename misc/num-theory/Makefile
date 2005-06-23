############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC)


EXES= angular buddy hyperg ndist series

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


angular:	angular.o $(LIB)/gcf.o $(LIB)/modular.o $(LIB)/moebius.o $(LIB)/totient.o
buddy:	buddy.o
hyperg:	hyperg.o $(LIB)/modular.o
ndist:	ndist.o $(LIB)/gcf.o $(LIB)/totient.o
series:	series.o divisor.o $(LIB)/gcf.o $(LIB)/moebius.o $(LIB)/totient.o

ndist.o:	ndist.c $(INC)/gcf.h $(INC)/totient.h 

divisor.o:	divisor.c divisor.h
