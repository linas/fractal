############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= fdist gappy sway

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) $(ENUMS)

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< $(LIB)/Farey.o $(LIB)/FareyTree.o -lm -lstdc++


gappy.o:  gappy.C $(LIB)/Farey.o $(INC)/Farey.h
sway.o:  sway.C $(LIB)/Farey.o $(INC)/Farey.h
fdist.o:  fdist.C $(LIB)/FareyTree.o $(INC)/FareyTree.h
