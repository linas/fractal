############### Symbol Definitions ###################### 

INCLUDES =	-I../tools/inc
LIBS =		/usr/lib

CC = cc -Wall -g -O2 $(INCLUDES)

EXES =  stalk stalkmov mstop migrate measure offset \
        lyapunov next squige age phase orig whack \
        circout circin circmov brat gapper \
        abs clamp takelog reclamp dump renorm \
        draw dx image movie paste stats sum \
        angle curl diff div card mandel morph \
        moire tong automatic wind winds ray \
        manvert inf decide cutoff

all:    ray automatic radius brat draw dx image


clean:
	rm -f a.out  *.o glop* junk* *.flo  *.txt *.gif  $(EXES)

.c.o:
	$(CC) -c $<

.C.o:
	$(CC) -c $<

.o:
	$(CC) -o $* $< -lm

image.o: image.h util.h
gap-tongue.o: image.h util.h
man.o: image.h util.h
util.o:	util.h

################## Dependencies and compile/link #########
#
gapper inf manvert stalk stalkmov mstop migrate measure offset lyapunov next squige age phase orig whack circout circin circmov brat: brat.o  util.o
	$(CC) -o brat brat.o util.o ../tools/lib/Farey.o ../tools/lib/FareyTree.o -lm -lstdc++
	ln -f brat decide
	ln -f brat manvert
	ln -f brat stalk
	ln -f brat mstop
	ln -f brat measure
	ln -f brat offset
	ln -f brat lyapunov
	ln -f brat next
	ln -f brat squige
	ln -f brat age
	ln -f brat migrate
	ln -f brat phase
	ln -f brat orig
	ln -f brat whack
	ln -f brat circout
	ln -f brat circin
	ln -f brat circmov
	ln -f brat stalkmov
	ln -f brat inf
	ln -f brat cutoff
	ln -f brat wind
	ln -f brat winds
	ln -f brat morph
	ln -f brat gapper
	ln -f brat tong

card: card.c  util.o
	$(CC) -o card card.c util.o -lm

draw: draw.c  util.o
	$(CC) -o draw -I /home/linas/linas/fractal-tools/inc draw.c util.o /home/linas/linas/fractal-tools/lib/Farey.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o /src/tools/lib/Farey.o -lm

sum diff paste angle curl div: sum.o  util.o opers.o
	$(CC) -o sum sum.o util.o opers.o -lm
	ln -f sum diff
	ln -f sum paste
	ln -f sum angle
	ln -f sum curl
	ln -f sum div

dx:	dx.c  util.o opers.o
	$(CC)  -o dx dx.c util.o opers.o -lm

abs clamp takelog reclamp dump renorm:	renorm.c  util.o opers.o
	$(CC)  -o renorm renorm.c opers.o util.o -lm
	ln -f renorm abs
	ln -f renorm clamp
	ln -f renorm takelog
	ln -f renorm reclamp
	ln -f renorm dump

stats:	stats.c  util.o opers.o
	$(CC)  -o stats stats.c opers.o util.o -lm

image: image.o  gap-tongue.o man.o util.o ../tools/lib/Farey.o
	$(CC) -o image $^ -lm -lstdc++
	ln -f image mandel
	ln -f image moire

movie: movie.o  man.o util.o
	$(CC) -o movie movie.o man.o util.o -lm

radius: radius.o 

automatic: radius
	ln -f radius automatic

ray:	ray.o
