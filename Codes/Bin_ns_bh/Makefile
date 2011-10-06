SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g
LIB = -L$(HOME_LORENE)/Lib -llorene -llorenef77

.SUFFIXES : .o .C

SRC1 = coal_ns_bh.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = init_ns_bh.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = lit_bin_ns_bh.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ1) $(LIB) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ2) $(LIB) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE3): $(OBJ3)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ3) $(LIB) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
					
.C.o:
	$(CXX)  -c $(CXXFLAGS) $(INC) $<

all:  $(EXE1) $(EXE2) $(EXE3) 

uninstall:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3)
	rm -f $(EXE1) $(EXE2) $(EXE3)

