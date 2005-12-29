
EXES= weyl

all: $(EXES)


weyl:	weyl.o

.o:
	cc -o $@ $< -lm

clean:
	rm *.o junk x

realclean: clean
	rm $(EXES)  *.png


