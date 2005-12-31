############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= c-takagi cantor-dist conj dyadic entropy exacto fdist feigen \
      fiddle fliptree gappy gkw ising mod moebius \
      mute orbit orbit-s padic \
      pair pegs quest rand-measure rdist summy swap sway sym \
      sym-fun takagi takagi-n tent-map theta tree vanish zeta

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf gap.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< $(LIB)/modular.o $(LIB)/gcf.o $(LIB)/Farey.o $(LIB)/FareyTree.o -lgsl -lcblas -lgmpxx -lgmp -lm -lstdc++ 


c-takagi.o:	c-takagi.c
takagi-n.o:	takagi-n.C
cantor-dist.o:	cantor-dist.C
conj.o:	conj.c
dyadic.o:	dyadic.C
entropy.o:  entropy.C $(LIB)/Farey.o $(INC)/Farey.h
exacto.o:  exacto.C $(LIB)/Farey.o $(INC)/Farey.h
fdist.o:  fdist.C $(LIB)/FareyTree.o $(INC)/FareyTree.h
feigen.o: feigen.C $(LIB)/Farey.o $(INC)/Farey.h
fiddle.o:  fiddle.C $(LIB)/Farey.o $(INC)/Farey.h
fliptree.o:  fliptree.C $(LIB)/Farey.o $(INC)/Farey.h
gappy.o:  gappy.C $(LIB)/Farey.o $(INC)/Farey.h
gkw.o:  gkw.C $(LIB)/Farey.o $(INC)/Farey.h
ising.o:  ising.C $(LIB)/Farey.o $(INC)/Farey.h
mod.o:  mod.C $(LIB)/modular.o $(INC)/modular.h
moebius.o:	moebius.C
mute.o:  mute.C $(LIB)/Farey.o $(INC)/Farey.h
orbit.o:	orbit.C
orbit-s.o:	orbit-s.C
padic.o:	padic.C
pegs.o:  pegs.C $(LIB)/gcf.o $(INC)/gcf.h
pair.o:  pair.C $(LIB)/Farey.o $(INC)/Farey.h
quest.o:  quest.C $(LIB)/Farey.o $(INC)/Farey.h
rand-measure.o:	rand-measure.C $(LIB)/gcf.o $(INC)/gcf.h
rdist.o:  rdist.C $(LIB)/Farey.o $(INC)/Farey.h
summy.o: summy.C $(LIB)/Farey.o $(INC)/Farey.h
swap.o:  swap.C $(LIB)/Farey.o $(INC)/Farey.h
sway.o:  sway.C $(LIB)/Farey.o $(INC)/Farey.h
sym.o: sym.C $(LIB)/Farey.o $(INC)/Farey.h
sym-fun.o: sym-fun.C $(LIB)/Farey.o $(INC)/Farey.h
takagi.o:	takagi.C
tent-map.o:	tent-map.C
theta.o:	theta.C
tree.o: tree.C $(LIB)/Farey.o $(INC)/Farey.h
vanish.o:	vanish.C
zeta.o:	zeta.C
