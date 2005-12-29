
EXES= ising takagi-gral

all: $(EXES)


ising:	ising.o
takagi-gral:	takagi-gral.o

.o:
	cc -o $@ $< -lm


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

