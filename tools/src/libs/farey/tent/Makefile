
SRCS= bernoulli.c logistic.c tent.c ftent.c

OBJS= $(SRCS:.c=.o)
EXES= $(SRCS:.c=)

INCLUDES= -I.. -I /src/tools/inc
LIBS= /src/tools/lib/gl.o /src/tools/lib/Farey.o

all: $(EXES)

# $(EXES):	$(OBJS)

$(OBJS):	$(SRCS)
	cc -c -O $(INCLUDES) $<

.o:
	cc -o $* $< $(LIBS) -lm -lX11


clean:
	rm -f *.o core tmp junk a.out

realclean: clean
	rm -f $(EXES)
