


EXES = stats classic wave

all: $(EXES)


classic:	classic.o ray.o
wave:	wave.o ray.o
stats: stats.o ray.o

.C.o:
	cc -g -O2 -c $<

.o:
	cc -o $* $^ -lm -lstdc++

clean:
	rm -f *.o

realclean: clean
	rm -f $(EXES)
