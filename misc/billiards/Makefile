


EXES = stats classic wave

all: $(EXES)


classic:	classic.o ray.o
	cc -o classic classic.o ray.o -lm


wave:	wave.o ray.o
	cc -o wave wave.o ray.o -lm

stats:	stats.o ray.o
	cc -o stats stats.o ray.o -lm

.C.o:
	g++ -g -O2 -c $<

clean:
	rm -f *.o

realclean: clean
	rm -f $(EXES)
