
# build Gnu MP code
#


EXES= bound mp_zeta zero-cross

all:  $(EXES)

bound:	bound.o
mp_zeta:	mp_zeta.o
zero-cross:	zero-cross.o

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< -lgmp -lm
