############### Symbol Definitions ###################### 
TOP=../../../..

INC = $(TOP)/inc
LIB = $(TOP)/lib

# CC = /bin/fcc -f2 -O -I $(INC)
# CC = /bin/cc -O2 -O -I $(INC) -DAIX315
# CC = /bin/cc -g -I $(INC) -DAIX315
CC = /usr/bin/cc -O2 -I $(INC) -DLINUX

RM=	/bin/rm -f

EXES= newmap cmap dmap imap zmap
all: 	$(EXES)


.c.o:
	cc -O2 -c $< -I..

.C.o:
	cc -O2 -c $< -I..

.o:
	cc -o $* $< $(LIB)/Farey.o -lm -lstdc++

reiter:	reiter.o
	cc -o reiter reiter.o  ../../../../lib/Farey.o -lm


################## Dependencies and compile/link #########
# newmap
# 	(Jan 1994) excercise Farey Number converter
#
prodmap: prodmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o prodmap prodmap.c $(LIB)/Farey.o -lm
	$(RM) prodmap.o

# cmap
# 	(Jan 1994) excercise Farey Number converter
#
cmap: cmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o cmap cmap.c $(LIB)/Farey.o -lm
	$(RM) cmap.o

# dmap
# 	(Jan 1994) excercise Farey Number converter
#
dmap: dmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o dmap dmap.c $(LIB)/Farey.o -lm
	$(RM) dmap.o

# imap
# 	(Jan 1994) excercise Farey Number converter
#
imap: imap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o imap imap.c $(LIB)/Farey.o -lm
	$(RM) imap.o

# flipmap
# 	(Jan 1994) excercise Farey Number converter
#
flipmap: flipmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o flipmap flipmap.c $(LIB)/Farey.o -lm
	$(RM) flipmap.o

symmap: symmap.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o symmap symmap.c $(LIB)/Farey.o -lm
	$(RM) symmap.o

clean:
	$(RM) a.out core *.o tmp junk glop

realclean:	clean
	$(RM) $(EXES) *.flo *.gif
	rm -f pixmaps


###########################################################

# ------------------------------------------------------------
enums:
	elist 5000 > enums

enums.ps:
	filter enums i f | xgraph -t "Exp'd Continued Fractions"

# ------------------------------------------------------------
fnums:
	flist 5000 > fnums

fnums.ps:
	filter fnums i f | xgraph -t "Farey Numbers"

# ------------------------------------------------------------

gnums:
	glist 3000 0.01 > gnum0.01 
	glist 1000 0.1 > gnum0.1 
	glist 1000 0.4 > gnum0.4 
	glist 1000 0.6 > gnum0.6 
	glist 1000 0.69314718 > gnum.ln2 
	glist 1000 0.8 > gnum0.8
	glist 1000 1.0 > gnum1.0

gnums.ps:
	filter gnum0.01 i f | xgraph -t "Farey Number with t=0.01" &
	filter gnum0.1 i f | xgraph -t "Farey Number with t=0.1" &
	filter gnum0.4 i f | xgraph -t "Farey Number with t=0.4" &
	filter gnum0.6 i f | xgraph -t "Farey Number with t=0.6" & 
	filter gnum.ln2 i f | xgraph -t "Farey Number with t=ln2" &
	filter gnum0.8 i f | xgraph -t "Farey Number with t=0.8" &
	filter gnum1.0 i f | xgraph -t "Farey Number with t=1.0" &

# ------------------------------------------------------------

znums:
#	zlist 3000 -1.2 > znum-1.2
#	zlist 3000 -1.0 > znum-1.0
#	zlist 3000 -0.8 > znum-0.8
#	zlist 3000 -0.5 > znum-0.5
	zlist 3000 -0.3 > znum-0.3
	zlist 3000 -0.25 > znum-0.25
#	zlist 3000 -0.2 > znum-0.2
#	zlist 3000 -0.03 > znum-0.03
#	zlist 3000 0.03 > znum0.03
#	zlist 3000 0.2 > znum0.2
#	zlist 3000 0.5 > znum0.5
#	zlist 3000 0.8 > znum0.8
#	zlist 3000 1.0 > znum1.0
#	zlist 3000 1.2 > znum1.2
#	zlist 3000 1.6 > znum1.6
#	zlist 3000 2.0 > znum2.0
#	zlist 3000 4.0 > znum4.0
#	zlist 3000 10.0 > znum10.0

znums.ps:
#	filter znum-1.2 i f | xgraph -t "Real Number with z= -1.2" &
#	filter znum-1.0 i f | xgraph -t "Real Number with z= -1.0" &
#	filter znum-0.8 i f | xgraph -t "Real Number with z= -0.8" &
#	filter znum-0.5 i f | xgraph -t "Real Number with z= -0.5" &
	filter znum-0.3 i f | xgraph -t "Real Number with z= -0.3" &
	filter znum-0.25 i f | xgraph -t "Real Number with z= -0.25" &
#	filter znum-0.2 i f | xgraph -t "Real Number with z= -0.2" &
#	filter znum-0.03 i f | xgraph -t "Real Number with z= -0.03" &
#	filter znum0.03 i f | xgraph -t "Real Number with z=0.03" &
#	filter znum0.2 i f | xgraph -t "Real Number with z=0.2" &
#	filter znum0.5 i f | xgraph -t "Real Number with z=0.5" &
#	filter znum0.8 i f | xgraph -t "Real Number with z=0.8" &
#	filter znum1.0 i f | xgraph -t "Real Number with z=1.0" &
#	filter znum1.2 i f | xgraph -t "Real Number with z=1.2" &
#	filter znum1.6 i f | xgraph -t "Real Number with z=1.6" &
#	filter znum2.0 i f | xgraph -t "Real Number with z=2.0" &
#	filter znum4.0 i f | xgraph -t "Real Number with z=4.0" &
#	filter znum10.0 i f | xgraph -t "Real Number with z=10.0" &

# ------------------------------------------------------------

zinvs:
	zinv 3000 0.03 > zinv0.03
	zinv 3000 0.2 > zinv0.2
	zinv 3000 0.5 > zinv0.5
	zinv 3000 2.0 > zinv2.0
	zinv 3000 4.0 > zinv4.0
	zinv 3000 10.0 > zinv10.0

zinvs.ps:
	filter zinv0.03 i f | xgraph -t "Reciprocal Real with z=0.03" &
	filter zinv0.2 i f | xgraph -t "Reciprocal Real with z=0.2" &
	filter zinv0.5 i f | xgraph -t "Reciprocal Real with z=0.5" &
	filter zinv2.0 i f | xgraph -t "Reciprocal Real with z=2.0" &
	filter zinv4.0 i f | xgraph -t "Reciprocal Real with z=4.0" &
	filter zinv10.0 i f | xgraph -t "Reciprocal Real with z=10.0" &

# ------------------------------------------------------------

station.ps:
	cat station.z | xgraph -t "Pure Gaussian Noise Stationarity" -x "loz10 (t2-t1)" -y "loz10 Var (X(t2)-X(t1))"
	cat station.b | xgraph -t "Stationarity of Brownian Motion" -x "loz10 (t2-t1)" -y "loz10 Var (X(t2) - X(t1))"

# ------------------------------------------------------------

smooth:
	spectrum zauss.z i z > smooth.z
	spectrum brown.b i z > smooth.b


smooth.ps:
	cat smooth.z | xgraph -t "Smooth Loz Power Spectrum of Gaussian Noise" -x "loz10 (omeza)" -y "loz10 Power"
	cat smooth.b | xgraph -t "Smooth Loz Power Spectrum of Brownian Motion" -x "loz10 (omeza)" -y "loz10 Power"

