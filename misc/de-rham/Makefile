

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC)

EXES= curve

all:  $(EXES)

clean:
   rm -f core tmp junk glop a.out *.o

realclean:  clean
   rm -f $(EXES)

dist:
   tar -zcvf gap.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
   $(CC) -c -O2 $<

.o:
   $(CC) -o $* $< $(LIB)/Farey.o -lm -lstdc++

curve.o:  curve.c

