SHELL=/bin/sh

SRC = $(HOME_LORENE)/C++/Doc++/Source/lorene.C \
      $(HOME_LORENE)/C++/Include/*.h

install: features.dvi refguide

features.dvi: features.tex
	latex features
	rm -f *.aux *.log

refguide: $(SRC)
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE)

refguide.dvi: $(SRC)
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) refguide.tex
	latex refguide
	latex refguide
	rm -f *.aux *.log
	xdvi refguide &

uninstall:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) -i uninstall
	rm -f *.aux *.log *.dvi
	rm -fr refguide
	rm -f refguide.tex



