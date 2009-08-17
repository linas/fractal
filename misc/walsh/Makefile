############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC)


EXES= mob-count walsh walsh-mink

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf walsh.tgz *.c *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -O2 -Wall -I$(INC) $<

.C.o:
	$(CC) -c -O2 -Wall -I$(INC) $<

.o:
	$(CC) -o $* $^ -L$(LIB) -lfunc -lm -lstdc++


mob-count:	mob-count.o
walsh:	walsh.o
walsh-mink:	walsh-mink.o


