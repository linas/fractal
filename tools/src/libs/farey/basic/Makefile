############### Symbol Definitions ###################### 

TOP = ../../../..
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= fgraph

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) $(ENUMS)

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< $(LIB)/Farey.o -lm -lstdc++


fgraph.o:  fgraph.C $(LIB)/Farey.o $(INC)/Farey.h
