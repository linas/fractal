
# Make the submission to arxiv xxxxxxxxxxxmath.NT/0609775
#
# \includegraphics{figures/crossover-def.eps.gz}
# bmeps -leps2 what.png what.eps
# maketex -fixeps *.eps
# maketex -gzip *.eps
# 
#
#
# 
#
IMGS= poly-perf-cold-cache.png \
      poly-perf-warm-cache.png \
      hurwitz-perf-cold.png \
      hurwitz-perf-warm.png \
      periodic-xsection.png \
      periodic-50.png \
      hurwitz-non-leading.png \
      hurwitz-extended.png \
      poly-phase.png \
      polylog-0-sheet.png \
      polylog-1.2-14.png \
      polylog-p1-sheet.png \
      polylog-m1-sheet.png \
      polylog,-0,-2-sheet.png \
      polylog,-0,-3-sheet.png \
      polylog,1,1-sheet.png \
      polylog,1,-2-sheet.png \
      polylog,1,-3-sheet.png

all:
	cp hurwitz.tex hurwitz.latex
	/usr/bin/latex hurwitz.latex
	/usr/bin/bibtex hurwitz.aux

dist:
	tar -zcvf hurwitz.tar.gz hurwitz.latex hurwitz.aux hurwitz.bbl $(IMGS)

clean:
	rm -f *.log *.blg *.dvi
