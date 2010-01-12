
CC=gcc
CFLAGS= -g -O2 -Wall

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = alpha asubn bern bern_h binom-graph character \
	diagon diverge eicheck frobenius gsl-diag \
	harmonic rational-sum riemann stirling sum swap tryall valuate \
	zeta-sum zeta-takagi zeta-zero zero zyx

all: $(EXES)

# library files
ache.o: ache.h 

# executables
alpha.o:	ache.h $(INC)/binomial.h
asubn.o:	ache.h 
bern.o: 
bern_h.o: 
beta.o: ache.h zetafn.h
binom-graph.o: zetafn.h
diagon.o: zetafn.h
eicheck.o: ache.h zetafn.h
frobenius.o: zetafn.h
gsl-diag.o:	zetafn.h
matrix.o: matrix.h
riemann.o:	zetafn.h
striling.o:	$(INC)/binomial.h
sum.o:	ache.h zetafn.h
tryall.o:	zetafn.h
zeta-sum.o: zetafn.h
zeta-zero.o: zetafn.h
zyx.o:	zetafn.h

alpha: alpha.o ache.o
asubn: asubn.o ache.o
bern: bern.o
bern_h: bern_h.o
beta: beta.o ache.o 
binom-graph: binom-graph.o
character: character.o 
diagon: diagon.o matrix.o
diverge: diverge.o
eicheck: eicheck.o ache.o
frobenius: frobenius.o
gsl-diag: gsl-diag.o
harmonic: harmonic.o
rational-sum: rational-sum.o
riemann: riemann.o
stirling: stirling.o ache.o
sum: sum.o ache.o
swap: swap.o
tryall: tryall.o ache.o
valuate: valuate.o
zeta-sum: zeta-sum.o
zeta-takagi: zeta-takagi.o
zeta-zero: zeta-zero.o
zero: zero.o
zyx: zyx.o


.o:
	$(CC) -o $@ $^ -L$(LIB) -lfunc -lgsl -lgslcblas -lm

.c.o:
	$(CC) -c $(CFLAGS) -I $(INC) $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot

realclean: clean

dist:
	tar -zcvf x *.c *.h *.gplot README Makefile .cvsignore */*.c */*.h */*.gplot */Makefile */.cvsignore
