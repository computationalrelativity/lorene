SHELL=/bin/sh

SRC = lorene.C $(HOME_LORENE)/C++/Include/*.h 

install: features.dvi refguide

features.dvi: features.tex
	latex features
	rm -f *.aux *.log

refguide: $(SRC) doxyfile
	mv refguide html
	doxygen doxyfile
	mv html refguide

doxyfile: lorene.dox
	rm -f doxyfile
	cp lorene.dox doxyfile
	doxygen -u doxyfile

uninstall:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) -i uninstall
	rm -f *.aux *.log *.dvi doxyfile*
	rm -fr refguide
	rm -f refguide.tex



