

# Make the submission to arxiv math.NT/0609775

all:
	mv plouffe-ram.tex plouffe-ram.latex
	/usr/bin/latex plouffe-ram.latex
	/usr/bin/bibtex plouffe-ram.aux


dist:
	tar -zcvf plouffe-ram.tar.gz plouffe-ram.latex plouffe-ram.aux plouffe-ram.bbl *.eps
