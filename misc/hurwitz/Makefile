
EXES= hurwitz

all: $(EXES)

hurwitz:	hurwtiz.o

.c.o:
	cc -c -g -O2 $<

.o:
	cc -g -o $@ $^ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

