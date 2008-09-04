############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= f-wavelet fanal fdist fourier ftag janal jdist rdist

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf fdist.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -g -O2 $<

.o:
	$(CC) -o $* $< -L$(LIB) -lfunc -lgsl -lcblas -lgmpxx -lgmp -lm -lstdc++ 


f-wavelet.o:  f-wavelet.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fanal.o:  fanal.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fdist.o:  fdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fourier.o:  fourier.C $(LIB)/libfunc.a $(INC)/FareyTree.h
ftag.o:  ftag.C $(LIB)/libfunc.a
janal.o:  janal.C 
jdist.o:  jdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
rdist.o:  rdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
