
all: rexfer

rexfer.o: rexfer.c
rexfer: rexfer.o

.c.o:
	cc -c -O2 $<

.o:
	cc -o $* $<
