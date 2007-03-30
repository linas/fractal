
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
GRPHS= poly-perf-cold-cache.ps \
      poly-perf-warm-cache.ps \
      hurwitz-perf-cold.ps \
      hurwitz-perf-warm.ps \
      periodic-xsection.ps \
      monodromy.eps

IMGS= periodic-50.ps \
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

all: $(GRPHS) $(IMGS)
	cp hurwitz-arxiv.tex hurwitz.latex
	/usr/bin/latex hurwitz.latex
	/usr/bin/bibtex hurwitz.aux

arxiv:
	patch -p0 < arxiv.diff
	mv hurwitz.lyx hurwitz-arxiv.lyx
	cvs update
	lyx-qt hurwitz-arxiv.lyx -e latex
	mv hurwitz-arxiv.tex hurwitz.latex
	/usr/bin/latex hurwitz.latex
	/usr/bin/bibtex hurwitz.aux

dist: $(GRPHS) $(IMGS)
	tar -cvf hurwitz.tar hurwitz.latex hurwitz.aux hurwitz.bbl $(GRPHS) $(IMGS)
	ls -la hurwitz.tar
	gzip hurwitz.tar
	ls -la hurwitz.tar.gz

clean:
	rm -f *.log *.blg *.dvi

imgclean:
	rm -f $(IMGS)

.SUFFIXES: .ps .png

.png.ps:
	cat $^ |  pngtopnm | pnmscale 0.65 -nomix | pnmtopng -compression 9 > tmp.png
	convert tmp.png $*.ps
