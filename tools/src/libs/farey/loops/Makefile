############### Symbol Definitions ###################### 

INC = $(TOP)/inc
LIB = $(TOP)/lib

# CC = /bin/fcc -f2 -O -I $(INC)
CC = /bin/cc -O2 -O -I $(INC) -DAIX315
# CC = /bin/cc -g -I $(INC) -DAIX315

RM=	/bin/rm -f

EXES= loop
all: 	$(EXES)

################## Dependencies and compile/link #########
# loop
# 	(Jan 1994) excercise Farey Number converter
#
loop: loop.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o loop loop.c $(LIB)/Farey.o -lm -lgl
	$(RM) loop.o


clean:
	$(RM) a.out core *.o 

realclean:
	$(RM) a.out core *.o $(EXES)


###########################################################
