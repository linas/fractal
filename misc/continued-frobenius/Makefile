
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = alpha asubn bern bern_h binom-graph character \
	diagon diverge eicheck frobenius gsl-diag \
	harmonic riemann swap sum tryall valuate \
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
sum.o:	ache.h zetafn.h
riemann.o:	zetafn.h
tryall.o:	zetafn.h
zeta-sum.o: zetafn.h
zeta-zero.o: zetafn.h
zyx.o:	zetafn.h

alpha: alpha.o ache.o $(LIB)/binomial.o $(LIB)/harmonic.o
asubn: asubn.o ache.o $(LIB)/binomial.o $(LIB)/harmonic.o
bern: bern.o $(LIB)/binomial.o $(LIB)/harmonic.o $(LIB)/bernoulli.o
bern_h: bern_h.o $(LIB)/binomial.o $(LIB)/bernoulli.o
beta: beta.o ache.o 
binom-graph: binom-graph.o  $(LIB)/binomial.o
character: character.o 
diagon: diagon.o  $(LIB)/harmonic.o
diverge: diverge.o
eicheck: eicheck.o ache.o  $(LIB)/binomial.o $(LIB)/harmonic.o
frobenius: frobenius.o $(LIB)/binomial.o $(LIB)/harmonic.o
gsl-diag: gsl-diag.o $(LIB)/binomial.o $(LIB)/harmonic.o
harmonic: harmonic.o
riemann: riemann.o $(LIB)/binomial.o $(LIB)/harmonic.o
sum: sum.o ache.o  $(LIB)/binomial.o $(LIB)/harmonic.o
swap: swap.o
tryall: tryall.o ache.o $(LIB)/binomial.o $(LIB)/harmonic.o
valuate: valuate.o
zeta-sum: zeta-sum.o $(LIB)/binomial.o
zeta-takagi: zeta-takagi.o  $(LIB)/binomial.o $(LIB)/harmonic.o
zeta-zero: zeta-zero.o  $(LIB)/binomial.o $(LIB)/harmonic.o
zero: zero.o
zyx: zyx.o  $(LIB)/binomial.o $(LIB)/harmonic.o


.o:
	cc -o $@ $^ -lgsl -lgslcblas -lm

.c.o:
	cc -g -c -O3 -Wall -I $(INC) $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot

realclean: clean

dist:
	tar -zcvf x *.c *.h *.gplot README Makefile .cvsignore */*.c */*.h */*.gplot */Makefile */.cvsignore
