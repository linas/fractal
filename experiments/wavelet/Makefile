############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= question

all: 	$(EXES)

question.o:  question.C 

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

.C.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $< -lgsl -lcblas -lm -lstdc++ 


dist:
	tar -zcvf x *.C Makefile .cvsignore
