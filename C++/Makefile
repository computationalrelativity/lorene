SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		

BIN_TOOLS = $(HOME_LORENE)/bin
LIB	  = $(HOME_LORENE)/Lib

install: $(LIB) $(BIN_TOOLS) $(BIN_TOOLS)/pointCpointo
	cd Source/Base_val; $(MAKE)
	cd Source/Base_vect; $(MAKE)
	cd Source/Bhole; $(MAKE)
	cd Source/Bhole_binaire; $(MAKE)
	cd Source/Bin_ns_bh; $(MAKE)
	cd Source/Bin_ns_ncp; $(MAKE)
	cd Source/Binaire; $(MAKE)
	cd Source/Champ_cart; $(MAKE)
	cd Source/Change_var; $(MAKE)
	cd Source/Cmp; $(MAKE)
	cd Source/Connection; $(MAKE)
	cd Source/Coord; $(MAKE)
	cd Source/Eos; $(MAKE)
	cd Source/Etoile; $(MAKE)
	cd Source/Grille3d; $(MAKE)
	cd Source/Grille_cart; $(MAKE)
	cd Source/Itbl; $(MAKE)
	cd Source/Map; $(MAKE)
	cd Source/Matrice; $(MAKE)
	cd Source/Metric; $(MAKE)
	cd Source/Metrique; $(MAKE)
	cd Source/Mg3d; $(MAKE)
	cd Source/Mtbl; $(MAKE)
	cd Source/Mtbl_cf; $(MAKE)
	cd Source/Non_class_members/Coef; $(MAKE)
	cd Source/Non_class_members/Graphics; $(MAKE)
	cd Source/Non_class_members/Operators; $(MAKE)
	cd Source/Non_class_members/PDE; $(MAKE)
	cd Source/Non_class_members/Utilities; $(MAKE)
	cd Source/Ope_elementary; $(MAKE)
	cd Source/Param; $(MAKE)
	cd Source/Param_elliptic; $(MAKE)
	cd Source/QuadriD; $(MAKE)
	cd Source/Tbl; $(MAKE)
	cd Source/Tenseur; $(MAKE)
	cd Source/Tensor; $(MAKE)
	cd Source/Tensor/Scalar; $(MAKE)
	cd Source/Valencia; $(MAKE)
	cd Source/Valeur; $(MAKE)
	touch $(LIB)/date_mark

uninstall:
	rm -fr $(LIB)/liblorene.a
	rm -fr $(LIB)/liblorene_g.a
	rm -f $(BIN_TOOLS)/pointCpointo
	cd Source/Base_val; $(MAKE) -i uninstall
	cd Source/Base_vect; $(MAKE) -i uninstall
	cd Source/Bhole; $(MAKE) -i uninstall
	cd Source/Bhole_binaire; $(MAKE) -i uninstall
	cd Source/Bin_ns_bh; $(MAKE) -i uninstall
	cd Source/Bin_ns_ncp; $(MAKE) -i uninstall
	cd Source/Binaire; $(MAKE) -i uninstall
	cd Source/Champ_cart; $(MAKE) -i uninstall
	cd Source/Change_var; $(MAKE) -i uninstall
	cd Source/Cmp; $(MAKE) -i uninstall
	cd Source/Connection; $(MAKE) -i uninstall
	cd Source/Coord; $(MAKE) -i uninstall
	cd Source/Eos; $(MAKE) -i uninstall
	cd Source/Etoile; $(MAKE) -i uninstall
	cd Source/Grille3d; $(MAKE) -i uninstall
	cd Source/Grille_cart; $(MAKE) -i uninstall
	cd Source/Itbl; $(MAKE) -i uninstall
	cd Source/Map; $(MAKE) -i uninstall
	cd Source/Matrice; $(MAKE) -i uninstall
	cd Source/Metric; $(MAKE) -i uninstall
	cd Source/Metrique; $(MAKE) -i uninstall
	cd Source/Mg3d; $(MAKE) -i uninstall
	cd Source/Mtbl; $(MAKE) -i uninstall
	cd Source/Mtbl_cf; $(MAKE) -i uninstall
	cd Source/Non_class_members/Coef; $(MAKE) -i uninstall
	cd Source/Non_class_members/Graphics; $(MAKE) -i uninstall
	cd Source/Non_class_members/Operators; $(MAKE) -i uninstall
	cd Source/Non_class_members/PDE; $(MAKE) -i uninstall
	cd Source/Non_class_members/Utilities; $(MAKE) -i uninstall
	cd Source/Ope_elementary; $(MAKE) -i uninstall
	cd Source/Param; $(MAKE) -i uninstall
	cd Source/Param_elliptic; $(MAKE) -i uninstall
	cd Source/QuadriD; $(MAKE) -i uninstall
	cd Source/Tbl; $(MAKE) -i uninstall
	cd Source/Tenseur; $(MAKE) -i uninstall
	cd Source/Tensor; $(MAKE)  -i uninstall
	cd Source/Tensor/Scalar; $(MAKE)  -i uninstall
	cd Source/Valencia; $(MAKE) -i uninstall
	cd Source/Valeur; $(MAKE) -i uninstall

$(LIB):
	mkdir $@

$(BIN_TOOLS):	
	mkdir $@

$(BIN_TOOLS)/pointCpointo:
	$(CXX) $(CXXFLAGS) -o $@ $(HOME_LORENE)/Inst_tools/pointCpointo.C

