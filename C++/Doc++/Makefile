SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		

BIN_TOOLS = $(HOME_LORENE)/bin

install: $(BIN_TOOLS) $(BIN_TOOLS)/reorg_comments html

html: Source/lorene.C $(HOME_LORENE)/C++/Include/*.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) . 
	doc++ -p -d html lorene.C_r

tex: lorene.dvi

lorene.dvi: Source/lorene.C $(HOME_LORENE)/C++/Include/*.h
	$(BIN_TOOLS)/reorg_comments $?
	mv $(?:=_r) . 
	doc++ -p -t -o lorene.tex lorene.C_r
	latex lorene
	latex lorene
	xdvi lorene &

$(BIN_TOOLS):	
	mkdir $@

$(BIN_TOOLS)/reorg_comments:
	$(CXX) $(CXXFLAGS) -o $@ $(HOME_LORENE)/Inst_tools/reorg_comments.C

uninstall:
	rm -f $(BIN_TOOLS)/reorg_comments
	rm -f *.h_r
	rm -f *.C_r
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -fr html



	
