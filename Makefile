SHELL=/bin/sh

install: 
	cd Doc; $(MAKE) -i
	cd C++; $(MAKE)
	cd F77; $(MAKE)
	cd Export/C++; $(MAKE)

uninstall:
	cd Doc; $(MAKE) -i uninstall
	cd C++; $(MAKE) -i uninstall
	cd F77; $(MAKE) -i uninstall
	cd Export/C++; $(MAKE) -i uninstall
	rm -fr bin
	rm -fr Lib
	rm -f Test/*.o
	rm -f Test/test_fft Test/test_lapack Test/test_pgplot

test: 
	cd Test; $(MAKE)


