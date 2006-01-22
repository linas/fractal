

EXES= galois

all: $(EXES)


galois: galois.o

.o:
	cc -o $@ $< -lm

clean:
	rm -f *.o junk a.out x

realclean: clean
	rm -f $(EXES)  *.png

