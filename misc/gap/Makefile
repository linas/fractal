############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= entropy exacto fdist feigen fiddle fliptree gappy mute padic \
      rdist summy swap sway sym \
      sym-fun takagi theta tree vanish

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< $(LIB)/Farey.o $(LIB)/FareyTree.o -lgmpxx -lgmp -lm -lstdc++ 


entropy.o:  entropy.C $(LIB)/Farey.o $(INC)/Farey.h
exacto.o:  exacto.C $(LIB)/Farey.o $(INC)/Farey.h
fdist.o:  fdist.C $(LIB)/FareyTree.o $(INC)/FareyTree.h
feigen.o: feigen.C $(LIB)/Farey.o $(INC)/Farey.h
fiddle.o:  fiddle.C $(LIB)/Farey.o $(INC)/Farey.h
fliptree.o:  fliptree.C $(LIB)/Farey.o $(INC)/Farey.h
gappy.o:  gappy.C $(LIB)/Farey.o $(INC)/Farey.h
mute.o:  mute.C $(LIB)/Farey.o $(INC)/Farey.h
padic.o:	padic.C
sway.o:  sway.C $(LIB)/Farey.o $(INC)/Farey.h
swap.o:  swap.C $(LIB)/Farey.o $(INC)/Farey.h
rdist.o:  rdist.C $(LIB)/Farey.o $(INC)/Farey.h
summy.o: summy.C $(LIB)/Farey.o $(INC)/Farey.h
sym.o: sym.C $(LIB)/Farey.o $(INC)/Farey.h
sym-fun.o: sym-fun.C $(LIB)/Farey.o $(INC)/Farey.h
takagi.o:	takagi.C
theta.o:	theta.C
tree.o: tree.C $(LIB)/Farey.o $(INC)/Farey.h
vanish.o:	vanish.C
