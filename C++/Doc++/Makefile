SHELL=/bin/sh

include $(HOME_LORENE)/local_settings

BIN_TOOLS = $(HOME_LORENE)/bin

install: $(BIN_TOOLS) $(BIN_TOOLS)/reorg_comments $(HOME_LORENE)/Doc/refguide

$(HOME_LORENE)/Doc/refguide: Source/lorene.C $(HOME_LORENE)/C++/Include/*.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) .
	rm -fr $@
	doc++ -p -w -W -d $@ lorene.C_r

refguide.tex: Source/lorene.C $(HOME_LORENE)/C++/Include/*.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) .
	doc++ -p -w -W -t -o $(HOME_LORENE)/Doc/refguide.tex lorene.C_r

$(BIN_TOOLS):
	mkdir $@

$(BIN_TOOLS)/reorg_comments:
	$(CXX) $(CXXFLAGS) -o $@ $(HOME_LORENE)/Inst_tools/reorg_comments.C

uninstall:
	rm -f $(BIN_TOOLS)/reorg_comments
	rm -f *.h_r
	rm -f *.C_r
