TEX=tex
DVIPDF=dvipdfm

%.dvi:%.tex
	$(TEX) $<
	rm -f $(basename $<).log

%.pdf:%.dvi
	$(DVIPDF) $<

mpl.pdf: mpl.tex

clean:
	rm -f *.dvi *.log

veryclean: clean
	rm -f *.pdf
