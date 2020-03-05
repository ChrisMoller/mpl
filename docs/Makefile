TEX=tex
LATEX=latex
PDFLATEX=pdflatex
DVIPDF=dvips -f < $< | ps2pdf - $@


SOURCES = comparison-ops-table.tex \
          matrix-vector-ops-table.tex  \
          trigonometric-ops-table.tex \
          general-ops-table.tex \
          mpl.tex

%.pdf:%.ltx
	$(PDFLATEX) $<
	$(PDFLATEX) $<

%.dvi:%.ltx
	$(LATEX) $<
	rm -f $(basename $<).log

#%.dvi:%.tex
#	$(TEX) $<

%.pdf:%.dvi
	$(DVIPDF)

mpl.pdf: mpl.ltx

#mpl.pdf: mpl.dvi
#mpl.dvi: $(SOURCES)

clean:
	rm -f *.dvi *.log mpl.aux

veryclean: clean
	rm -f *.pdf
