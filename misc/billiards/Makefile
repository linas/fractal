


all: classic wave


classic:	classic.o ray.o
	cc -o classic classic.o ray.o -lm


wave:	wave.o ray.o
	cc -o wave wave.o ray.o -lm


.C.o:
	g++ -g -O2 -c $<
