############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

ANANT=../anant/src

CC = /usr/bin/cc -g -O2 -I$(INC) -I$(ANANT)

EXES= exact-bignum exact-series newton-series self swap yarh-lofi

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf yarh.tgz *.c *.C *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c $<

.C.o:
	$(CC) -c $<

.o:
	$(CC) -o $* $< -L$(LIB) -L$(ANANT) -lanant -lfunc -lgmp -ldb -lstdc++ -lm

exact-bignum: exact-bignum.o
exact-series: exact-series.o
newton-series: newton-series.o
self:  self.o
swap:  swap.o
yarh-lofi: yarh-lofi.o

swap.o:  swap.C $(INC)/Farey.h
