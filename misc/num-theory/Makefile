############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -Wall -I $(INC)


EXES= angular buddy divisor fourier hyperg interval liouville \
	mertens ndist series von_mangoldt

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf num.tgz *.c *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -O2 -I$(INC) $<

.C.o:
	$(CC) -c -O2 -I$(INC) $<

.o:
	$(CC) -o $* $^ -L$(LIB) -lfunc -lm -lstdc++


angular:	angular.o 
buddy:	buddy.o
divisor:	divisor.o
fourier:	fourier.o
hyperg:	hyperg.o
interval:	interval.o
liouville:	liouville.o
mertens:	mertens.o 
ndist:	ndist.o
series:	series.o  
von_mangoldt:	von_mangoldt.o

ndist.o:	ndist.c $(INC)/gcf.h $(INC)/totient.h 

