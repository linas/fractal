
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
IMGS= poly-perf-cold-cache.ps \
      poly-perf-warm-cache.ps \
      hurwitz-perf-cold.ps \
      hurwitz-perf-warm.ps \
      periodic-xsection.ps \
      periodic-50.ps \
      hurwitz-non-leading.ps \
      hurwitz-extended.ps \
      poly-phase.ps \
      polylog-0-sheet.ps \
      polylog-1.2-14.ps \
      polylog-p1-sheet.ps \
      polylog-m1-sheet.ps \
      polylog,-0,-2-sheet.ps \
      polylog,-0,-3-sheet.ps \
      polylog,1,1-sheet.ps

all:
	cp hurwitz-arxiv.tex hurwitz.latex
	/usr/bin/latex hurwitz.latex
	/usr/bin/bibtex hurwitz.aux

dist:
	tar -zcvf hurwitz.tar.gz hurwitz.latex hurwitz.aux hurwitz.bbl $(IMGS)
	ls -la hurwitz.tar.gz

clean:
	rm -f *.log *.blg *.dvi
