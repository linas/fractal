
all: rexfer

rexfer.o: rexfer.c
rexfer: rexfer.o

.c.o:
	cc -c -std=c99 -O2 $<

.o:
	cc -o $* $< -lm
