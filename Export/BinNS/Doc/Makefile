SHELL=/bin/sh

include $(HOME_LORENE)/local_settings

BIN_TOOLS = $(HOME_LORENE)/bin

INCEXP    = $(HOME_LORENE)/Export/C++/Include

install: $(BIN_TOOLS) $(BIN_TOOLS)/reorg_comments html documentation.tex

html: Source/documentation.C $(INCEXP)/bin_ns.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) .
	rm -fr $@
	doc++ -p -w -W -d $@ documentation.C_r

documentation.tex: Source/documentation.C $(INCEXP)/bin_ns.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) .
	doc++ -p -w -W -t -o $@ documentation.C_r


$(BIN_TOOLS):
	mkdir $@

$(BIN_TOOLS)/reorg_comments:
	$(CXX) $(CXXFLAGS) -o $@ $(HOME_LORENE)/Inst_tools/reorg_comments.C

uninstall:
	rm -f *.h_r
	rm -f *.C_r
	rm -f *.tex
	rm -fr html

