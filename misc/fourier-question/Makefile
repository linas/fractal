
EXES= gral-matrix gral-simple

all: $(EXES)


gral-matrix: gral-matrix.o
gral-simple: gral-simple.o

.c.o:
	cc -c -g -O2 $<

.o:
	cc -g -o $@ $^ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

