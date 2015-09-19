
all: fusc stern

fusc: fusc.o
stern: stern.o

.o:
	cc -o $@ $^ -lm -lstdc++


.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
