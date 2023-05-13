
EXES= weyl

all: $(EXES)


weyl:	weyl.o

.o:
	cc -o $@ $< -lm

clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png


