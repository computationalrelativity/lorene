SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

SRC1 = coal_bh.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = lit_bh.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = init_bh.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

SRC4 = prepare_seq.C
OBJ4 = $(SRC4:.C=.o)
EXE4 = $(SRC4:.C=)

SRC5 = init_coal.C
OBJ5 = $(SRC5:.C=.o)
EXE5 = $(SRC5:.C=)

SRC6 = coal_seq.C
OBJ6 = $(SRC6:.C=.o)
EXE6 = $(SRC6:.C=)

SRC7 = coal_highres.C
OBJ7 = $(SRC7:.C=.o)
EXE7 = $(SRC7:.C=)

$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE3): $(OBJ3)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE4): $(OBJ4)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE5): $(OBJ5)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(LIB_G) $(LIB_LAPACK) \
                                        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE6): $(OBJ6)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ6) $(LIB_G) $(LIB_LAPACK) \
                                        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE7): $(OBJ7)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ7) $(LIB_G) $(LIB_LAPACK) \
                                        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7)
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7)

