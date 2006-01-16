SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a

.SUFFIXES : .o .C

EXE1 = isolhor
SRC1 = isolhor.C
OBJ1 = $(SRC1:.C=.o)

EXE2 = lit_ih
SRC2 = lit_ih.C
OBJ2 = $(SRC2:.C=.o)

EXE3 = prepare_seq
SRC3 = prepare_seq.C
OBJ3 = $(SRC3:.C=.o)

EXE4 = isolhor_spher
SRC4 = isolhor_spher.C
OBJ4 = $(SRC4:.C=.o)


$(EXE1): $(OBJ1) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)
$(EXE2): $(OBJ2) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)
$(EXE3): $(OBJ3) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)
$(EXE4): $(OBJ4) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)





.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4)

