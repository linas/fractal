
EXES = bern bern_h diagon diverge eicheck frobenius gsl-diag \
	riemann sum tryall valuate zeta-sum zeta-zero zero

all: $(EXES)

# library files
ache.o: ache.h zetafn.h
zetafn.o: zetafn.h

# executables
bern.o: zetafn.h
bern_h.o: zetafn.h
beta.o: ache.h zetafn.h
diagon.o: zetafn.h
eicheck.o: ache.h zetafn.h
frobenius.o: zetafn.h
gsl-diag.o:	zetafn.h
sum.o:	ache.h zetafn.h
riemann.o:	zetafn.h
tryall.o:	zetafn.h
zeta-sum.o: zetafn.h
zeta-zero.o: zetafn.h

bern: bern.o zetafn.o
bern_h: bern_h.o zetafn.o
beta: beta.o ache.o zetafn.o
diagon: diagon.o zetafn.o
diverge: diverge.o
eicheck: eicheck.o ache.o zetafn.o 
frobenius: frobenius.o zetafn.o
gsl-diag: gsl-diag.o zetafn.o 
riemann: riemann.o zetafn.o
sum: sum.o ache.o zetafn.o
tryall: tryall.o ache.o zetafn.o
valuate: valuate.o
zeta-sum: zeta-sum.o zetafn.o
zeta-zero: zeta-zero.o zetafn.o 
zero: zero.o


.o:
	cc -o $@ $^ -lgsl -lgslcblas -lm

.c.o:
	cc -g -c -O3 -Wall $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot
