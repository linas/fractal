


EXES = stats classic wave

all: $(EXES)


classic:	classic.o ray.o
	cc -o classic classic.o ray.o -lm -lstdc++


wave:	wave.o ray.o
	cc -o wave wave.o ray.o -lm -lstdc++

stats:	stats.o ray.o
	cc -o stats stats.o ray.o -lm -lstdc++

.C.o:
	cc -g -O2 -c $<

clean:
	rm -f *.o

realclean: clean
	rm -f $(EXES)
