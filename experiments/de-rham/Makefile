

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC)

EXES= affine curve de-rham

all:  $(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES)

dist:
	tar -zcvf gap.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -O2 $<

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< -L$(LIB) -lfunc -lm -lstdc++

affine.o:  affine.c
curve.o:  curve.C
de-rham.o:  de-rham.c

