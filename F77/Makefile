SHELL=/bin/sh

include $(HOME_LORENE)/local_settings	

LIB	  = $(HOME_LORENE)/Lib

install: $(LIB) $(LIB)/liblorenef77_g.a $(LIB)/liblorenef77.a

$(LIB)/liblorenef77_g.a: Source/*.f Source/Poisson2d/*.f
	cd Source; $(MAKE) $@
	cd Source/Poisson2d; $(MAKE) $@
	ar -rv $@ $(?:.f=.o)
	$(RANLIB) $@
	rm -f $(?:.f=.o)

$(LIB)/liblorenef77.a: Source/*.f Source/Poisson2d/*.f
	cd Source; $(MAKE) $@
	cd Source/Poisson2d; $(MAKE) $@
	ar -rv $@ $(?:.f=.o)
	$(RANLIB) $@
	rm -f $(?:.f=.o)

uninstall:
	rm -f $(LIB)/liblorenef77.a
	rm -f $(LIB)/liblorenef77_g.a

$(LIB):
	mkdir $@
