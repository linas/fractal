############### Symbol Definitions ###################### 

INC = $(TOP)/inc
LIB = $(TOP)/lib

# CC = /bin/fcc -f2 -O -I $(INC)
CC = /bin/cc -O2 -O -I $(INC) -DAIX315
# CC = /bin/cc -g -I $(INC) -DAIX315

RM=	/bin/rm -f

EXES= cprod  inv
all: 	$(EXES)


.c.o:
	cc -O -c $< -I..

.C.o:
	cc -O -c $< -I..

cprod:	cprod.o
	cc -o cprod cprod.o  ../../../../lib/Farey.o -lm
	ln -f cprod cden

inv:	inv.o
	cc -o inv inv.o  ../../../../lib/Farey.o -lm

interval:	interval.o
	cc -o interval interval.o  ../../../../lib/Farey.o -lm

drip:	drip.o
	cc -o drip drip.o  ../../../../lib/Farey.o -lm

square:	square.o
	cc -o square square.o  ../../../../lib/Farey.o -lm

################## Dependencies and compile/link #########
# cmap
# 	(Jan 1994) excercise Farey Number converter
#

clean:
	$(RM) a.out core *.o tmp junk glop

realclean:	clean
	$(RM) $(EXES) *.flo *.gif
	cd pixmaps; rm -f *


###########################################################

stuff:
	cden xy-2.5.flo 500 500 -2.5
	cden xy-1.5.flo 500 500 -1.5
	cden xy-1.0.flo 500 500 -1.0
	cden xy-0.4.flo 500 500 -0.4
	cden xy+0.0.flo 500 500 +0.0
	cden xy+0.6.flo 500 500 0.6
	cden xy+1.1.flo 500 500 1.1
	cden xy+2.0.flo 500 500 2.0
	cden xy+2.8.flo 500 500 2.8
