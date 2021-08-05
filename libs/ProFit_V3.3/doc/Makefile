DEST = /acrm/www/html/software/profit/
all : ProFit.ps ProFit.pdf .web


ProFit.ps : ProFit.dvi
	dvips -o $@ $<

ProFit.dvi : ProFit.tex
	latex ProFit
	latex ProFit

ProFit.pdf : ProFit.tex ProFit.dvi
	pdflatex ProFit

.web : ProFit.pdf
	latex2html -local_icons -split 4 ProFit
	touch .web

install : 
	\cp -R ProFit/* $(DEST)/doc
	\cp ProFit.css $(DEST)/doc
#	\cp ProFit.pdf $(DEST)

clean : 
	\rm -rf ProFit.aux ProFit.dvi ProFit.log ProFit.ps ProFit
