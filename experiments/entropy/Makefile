
all: defacto entropy gauss-kuzmin

defacto: defacto.o
entropy: entropy.o
gauss-kuzmin: gauss-kuzmin.o

.o:
	g++ -o $@ $^ -L ../../tools/lib -lfunc -lm -lstdc++

.c.o:
	cc -O2 -c $<

.C.o:
	g++ -O2 -I ../../tools/inc -c $<

