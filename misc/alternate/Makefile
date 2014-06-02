
all: alt-binary alt-frac

alt-binary: alt-binary.o
	cc -o $@ $^ -lm

alt-frac: alt-frac.o
	cc -o $@ $^ -lm

.C.o:
	cc -c $^
