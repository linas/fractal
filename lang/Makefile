#
# FILE:
# Makefile


LIBSRCS=			\
	multihash.C		\
	wordhash.C

EXESRCS=			\
	analyze.C		\
	strip.C

LIBOBJS=$(LIBSRCS:.C=.o)
EXEOBJS=$(EXESRCS:.C=.o)
EXES=$(EXESRCS:.C=)

.C.o:
	cc -c -g -O2 $< 

all: $(EXES)

analyze: analyze.o $(LIBOBJS)
	cc analyze.o -o analyze $(LIBOBJS) -lm

strip: strip.o 
	cc strip.o -o strip 

clean:
	rm -f junk glop tmp a.out core
	
realclean: clean
	rm -f $(LIBOBJS) $(EXEOBJS) $(EXES)
	makedepend
	rm -f Makefile.bak

depend:
	makedepend -- $(CLAGS) $(DEFINES) $(INCLUDES) \
	-I /usr/lib/g++-include			\
	 -I /usr/lib/gcc-lib/i486-linux/2.5.8/include -- \
	$(LIBSRCS) $(EXESRCS)

# ------------------------------------------------------------
# DO NOT DELETE THIS LINE -- makedepend depends on it.

multihash.o: multihash.h pairhash.h /usr/include/typedefs.h
multihash.o: /usr/include/stdlib.h /usr/include/features.h
multihash.o: /usr/include/sys/cdefs.h /usr/include/stddef.h
multihash.o: /usr/include/errno.h /usr/include/linux/errno.h
multihash.o: /usr/lib/gcc-lib/i486-linux/2.5.8/include/float.h
multihash.o: /usr/include/signal.h /usr/include/sys/types.h
multihash.o: /usr/include/linux/types.h /usr/include/asm/types.h
multihash.o: /usr/include/linux/signal.h /usr/include/sys/time.h
multihash.o: /usr/include/linux/time.h /usr/include/sys/time.h config.h
multihash.o: concord.h pairhash.C /usr/include/math.h /usr/include/huge_val.h
multihash.o: /usr/include/endian.h /usr/include/bytesex.h
multihash.o: /usr/include/values.h /usr/include/stdio.h /usr/include/libio.h
multihash.o: /usr/include/_G_config.h top.h concord.C
wordhash.o: /usr/include/stdio.h /usr/include/libio.h
wordhash.o: /usr/include/_G_config.h top.h wordhash.h config.h
analyze.o: /usr/include/ctype.h /usr/include/features.h
analyze.o: /usr/include/sys/cdefs.h /usr/include/stdio.h /usr/include/libio.h
analyze.o: /usr/include/_G_config.h /usr/include/string.h
analyze.o: /usr/include/stddef.h multihash.h pairhash.h
analyze.o: /usr/include/typedefs.h /usr/include/stdlib.h /usr/include/errno.h
analyze.o: /usr/include/linux/errno.h
analyze.o: /usr/lib/gcc-lib/i486-linux/2.5.8/include/float.h
analyze.o: /usr/include/signal.h /usr/include/sys/types.h
analyze.o: /usr/include/linux/types.h /usr/include/asm/types.h
analyze.o: /usr/include/linux/signal.h /usr/include/sys/time.h
analyze.o: /usr/include/linux/time.h /usr/include/sys/time.h config.h
analyze.o: concord.h wordhash.h
strip.o: /usr/include/ctype.h /usr/include/features.h
strip.o: /usr/include/sys/cdefs.h /usr/include/stdio.h /usr/include/libio.h
strip.o: /usr/include/_G_config.h /usr/include/string.h /usr/include/stddef.h
