
all: gauss-kuzmin

gauss-kuzmin: gauss-kuzmin.o

.o:
	cc -o $@ $^ -lm

.c.o:
	cc -O2 -c $<

