############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= swap yarh-lofi

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf yarh.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< -L$(LIB) -lfunc -lstdc++ -lm

swap:  swap.o $(LIB)/libfunc.a

swap.o:  swap.C $(INC)/Farey.h
yarh-lofi: yarh-lofi.o
