SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		

BIN_TOOLS = $(HOME_LORENE)/bin
LIB	  = $(HOME_LORENE)/Lib

install: $(LIB) $(BIN_TOOLS) $(BIN_TOOLS)/pointCpointo
	cd Source; $(MAKE)

uninstall:
	rm -fr $(LIB)/liblorene_export.a
	rm -fr $(LIB)/liblorene_export_g.a
	cd Source; $(MAKE) -i uninstall

$(LIB):
	mkdir $@
	touch date_mark

$(BIN_TOOLS):	
	mkdir $@

$(BIN_TOOLS)/pointCpointo:
	$(CXX) $(CXXFLAGS) -o $@ $(HOME_LORENE)/Inst_tools/pointCpointo.C

