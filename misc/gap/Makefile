############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= exacto fdist feigen gappy padic rdist swap sway sym sym-fun 

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) $(ENUMS)

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< $(LIB)/Farey.o $(LIB)/FareyTree.o -lm -lstdc++


exacto.o:  exacto.C $(LIB)/Farey.o $(INC)/Farey.h
fdist.o:  fdist.C $(LIB)/FareyTree.o $(INC)/FareyTree.h
gappy.o:  gappy.C $(LIB)/Farey.o $(INC)/Farey.h
padic.o:	padic.C
sway.o:  sway.C $(LIB)/Farey.o $(INC)/Farey.h
swap.o:  swap.C $(LIB)/Farey.o $(INC)/Farey.h
rdist.o:  rdist.C $(LIB)/Farey.o $(INC)/Farey.h
feigen.o: feigen.C $(LIB)/Farey.o $(INC)/Farey.h
sym.o: sym.C $(LIB)/Farey.o $(INC)/Farey.h
sym-fun.o: sym-fun.C $(LIB)/Farey.o $(INC)/Farey.h
