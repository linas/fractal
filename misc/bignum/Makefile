
# build Gnu MP code
#


EXES= bound bs bsubn cache-fill cont_frac dsubn minima mp_zeta_test \
      stieltjes trunc zero-cross

all:  $(EXES)

bound:	bound.o
cache-fill:	cache-fill.o mp_zeta.o db-cache.o
cont_frac:	cont_frac.o
dsubn:	dsubn.o mp_zeta.o db-cache.o
mp_zeta_test:	mp_zeta_test.o mp_zeta.o db-cache.o
stieltjes:	stieltjes.o mp_zeta.o db-cache.o
trunc:	trunc.o
zero-cross:	zero-cross.o

bs.o: mp_zeta.h
bsubn.o: mp_zeta.h
db-cache.o: db-cache.h
dsubn.o: mp_zeta.h
mp_zeta.o: mp_zeta.h db-cache.h
mp_zeta_test.o: mp_zeta.h

bs:	bs.o mp_zeta.o db-cache.o
	$(CC) -o bs $^ -lgmp -lgsl -lcblas -ldb -lm -lc

bsubn:	bsubn.o mp_zeta.o db-cache.o
	$(CC) -o bsubn $^ -lgmp -lgsl -lcblas -ldb -lm -lc

minima: minima.o
	$(CC) -o minima $< -lgsl -lcblas -lm -lc

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.h *.gplot *.txt bzero.dat Makefile .cvsignore

.c.o:
	$(CC) -c -g -O2 -Wall $<

.o:
	$(CC) -g -o $* $^ -lgmp -ldb -lm
