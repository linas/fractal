
############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -std=gnu11 -Wall -I $(INC)

EXES= ffbasis

all:  $(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) 
dist:
	tar -zcvf num.tgz *.c *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -O2 -I$(INC) $<

.C.o:
	$(CC) -c -O2 -I$(INC) $<

.o:
	$(CC) -o $* $^ -L$(LIB) -lfunc -lm -lstdc++


falling.o: falling.h $(INC)/binomial.h $(INC)/cache.h
ffbasis.o: falling.h $(INC)/binomial.h $(INC)/cache.h
ffbasis: falling.o ffbasis.o


