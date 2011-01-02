

# Make the submission to arxiv math.NT/0609775

all:

plouffe-ram.bbl: plouffe-ram.tex
	/usr/bin/latex plouffe-ram.tex
	/usr/bin/bibtex plouffe-ram.aux

plouffe-ram.tex: plouffe-ram.lyx
	lyx plouffe-ram.lyx -e latex

dist: plouffe-ram.tex plouffe-ram.bbl
	tar -zcvf plouffe-ram.tar.gz plouffe-ram.tex plouffe-ram.bbl *.eps

clean:
	rm -f *.log *.blg *.dvi
