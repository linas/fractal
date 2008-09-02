############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= zeta-fourier

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf fredholm.tgz *.c *.C *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -g -O2 $<

.o:
	$(CC) -o $* $< -L$(LIB) -lfunc -lgsl -lcblas -lm -lstdc++ 


zeta-fourier.o:  zeta-fourier.c $(LIB)/libfunc.a $(INC)/FareyTree.h
