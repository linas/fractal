############### Symbol Definitions ###################### 

INCLUDES =	-I/usr/include/X11 
LIBS =		/usr/lib
# XLIB =		$(LIBS)/libX11.a \
#		$(LIBS)/libbsd.a \
#		$(LIBS)/libsock.a

XLIB =		$(LIBS)/libX11.a 

CC = cc -Wall -g -O2 $(INCLUDES)

EXES =  stalk stalkmov mstop migrate measure offset \
        lyapunov next squige age phase orig whack \
        circout circin circmov brat \
        abs clamp takelog reclamp dump renorm \
        draw dx image movie paste stats sum \
        angle curl diff div card mandel

all:    radius brat draw dx image


clean:
	rm -f a.out  *.o glop* junk* *.flo  *.txt *.gif  $(EXES)

.c.o:
	$(CC) -c $<

.o:
	$(CC) -o $* $< -lm

################## Dependencies and compile/link #########
#
inf manvert stalk stalkmov mstop migrate measure offset lyapunov next squige age phase orig whack circout circin circmov brat: brat.c  util.o
	$(CC) -o brat brat.c util.o -lm
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

card: card.c  util.o
	$(CC) -o card card.c util.o -lm

draw: draw.c  util.o
	$(CC) -o draw -I /home/linas/linas/fractal-tools/inc draw.c util.o /home/linas/linas/fractal-tools/lib/Farey.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o /src/tools/lib/Farey.o -lm

sum: sum.c  util.o opers.o
	$(CC) -o sum sum.c util.o opers.o -lm
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

image: image.o  man.o util.o
	$(CC) -o image image.o man.o util.o -lm
	ln -f image mandel
	ln -f image moire

movie: movie.o  man.o util.o
	$(CC) -o movie movie.o man.o util.o -lm

radius: radius.o 

