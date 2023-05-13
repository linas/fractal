############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -Wall -I$(INC)

EXES= baire mob-count walsh walsh-mink

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf walsh.tgz *.c *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -std=c99 $<

.C.o:
	$(CC) -c $<

.o:
	$(CC) -o $* $^ -L$(LIB) -lfunc -lm -lstdc++


baire: baire.o
mob-count:	mob-count.o
walsh:	walsh.o
walsh-mink:	walsh-mink.o
