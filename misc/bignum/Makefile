
# build Gnu MP code
#


EXES= bound bs bsubn cont_frac minima stieltjes trunc zero-cross

all:  $(EXES)

bound:	bound.o
bs:	bs.o mp_zeta.o
bsubn:	bsubn.o mp_zeta.o
cont_frac:	cont_frac.o
stieltjes:	stieltjes.o mp_zeta.o
trunc:	trunc.o
zero-cross:	zero-cross.o

mp_zeta.o: mp_zeta.h
bsubn.o: mp_zeta.h

minima: minima.o
	$(CC) -o minima $< -lgsl -lcblas -lm -lc

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.h *.gplot *.txt Makefile .cvsignore

.c.o:
	$(CC) -c -g -O2 $<

.o:
	$(CC) -g -o $* $^ -lgmp -lm
