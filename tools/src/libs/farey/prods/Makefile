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

cprod:	cprod.o
	cc -o cprod cprod.o  ../../../../lib/Farey.o -lm
	ln -f cprod cden

inv:	inv.o
	cc -o inv inv.o  ../../../../lib/Farey.o -lm

################## Dependencies and compile/link #########
# cmap
# 	(Jan 1994) excercise Farey Number converter
#
cmap: cmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o cmap cmap.c $(LIB)/Farey.o -lm
	$(RM) cmap.o


clean:
	$(RM) a.out core *.o tmp junk glop

realclean:	clean
	$(RM) $(EXES) *.flo *.gif
	cd pixmaps; rm -f *


###########################################################

