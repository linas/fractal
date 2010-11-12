

all: zeta-fractal

zeta-fractal.o: zeta-fractal.c

.c.o:
	cc -O2 -c $<

.o:
	cc -o $* $^ -lm
