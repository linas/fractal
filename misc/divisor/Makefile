


############### Symbol Definitions ###################### 


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -Wall -I $(INC)


EXES= moebius

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf *.cc *.gplot Makefile .gitignore

.c.o:
	$(CC) -c -O2 -I$(INC) $<

.C.o:
	$(CC) -c -O2 -I$(INC) $<

.o:
	$(CC) -o $* $^ -L$(LIB) -lfunc -lm -lstdc++


moebius: moebius.o

moebius.o: moebius.cc $(INC)/totient.h


