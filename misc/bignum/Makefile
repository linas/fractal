
# build Gnu MP code
#


EXES= bound minima mp_zeta trunc zero-cross

all:  $(EXES)

bound:	bound.o
mp_zeta:	mp_zeta.o
trunc:	trunc.o
zero-cross:	zero-cross.o

minima: minima.o
	$(CC) -o minima $< -lgsl -lcblas -lm -lc

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.gplot *.txt Makefile .cvsignore

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< -lgmp -lm
