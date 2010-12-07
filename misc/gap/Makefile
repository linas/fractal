############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= c-takagi cantor-dist conj dyadic entropy exacto feigen \
      fiddle fliptree gappy gkw ising measure mod moebius \
      mute orbit orbit-s padic \
      pair pegs quest rand-measure summy sway sym \
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
	$(CC) -o $* $< -L$(LIB) -lfunc -lgsl -lgmpxx -lgmp -lm -lstdc++
#	$(CC) -o $* $< -L$(LIB) -lfunc -lgsl -lcblas -lgmpxx -lgmp -lm -lstdc++


c-takagi.o:	c-takagi.c
takagi-n.o:	takagi-n.C
cantor-dist.o:	cantor-dist.C
conj.o:	conj.c
dyadic.o:	dyadic.C
entropy.o:  entropy.C $(LIB)/libfunc.a $(INC)/Farey.h
exacto.o:  exacto.C $(LIB)/libfunc.a $(INC)/Farey.h
feigen.o: feigen.C $(LIB)/libfunc.a $(INC)/Farey.h
fiddle.o:  fiddle.C $(LIB)/libfunc.a $(INC)/Farey.h
fliptree.o:  fliptree.C $(LIB)/libfunc.a $(INC)/Farey.h
gappy.o:  gappy.C $(LIB)/libfunc.a $(INC)/Farey.h
gkw.o:  gkw.C $(LIB)/libfunc.a $(INC)/Farey.h
ising.o:  ising.C $(LIB)/libfunc.a $(INC)/Farey.h
measure.o:	measure.c
mod.o:  mod.C $(LIB)/libfunc.a $(INC)/modular.h
moebius.o:	moebius.C
mute.o:  mute.C $(LIB)/libfunc.a $(INC)/Farey.h
orbit.o:	orbit.C
orbit-s.o:	orbit-s.C
padic.o:	padic.C
pegs.o:  pegs.C $(LIB)/libfunc.a $(INC)/gcf.h
pair.o:  pair.C $(LIB)/libfunc.a $(INC)/Farey.h
quest.o:  quest.C $(LIB)/libfunc.a $(INC)/Farey.h
rand-measure.o:	rand-measure.C $(LIB)/libfunc.a $(INC)/gcf.h
summy.o: summy.C $(LIB)/libfunc.a $(INC)/Farey.h
sway.o:  sway.C $(LIB)/libfunc.a $(INC)/Farey.h
sym.o: sym.C $(LIB)/libfunc.a $(INC)/Farey.h
sym-fun.o: sym-fun.C $(LIB)/libfunc.a $(INC)/Farey.h
takagi.o:	takagi.C
tent-map.o:	tent-map.C
theta.o:	theta.C
tree.o: tree.C $(LIB)/libfunc.a $(INC)/Farey.h
vanish.o:	vanish.C
zeta.o:	zeta.C
