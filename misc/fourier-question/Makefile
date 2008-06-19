
EXES= gral-matrix

all: $(EXES)


gral-matrix: gral-matrix.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

