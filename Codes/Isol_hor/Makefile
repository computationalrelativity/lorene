SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a

.SUFFIXES : .o .C

EXE1 = isolhor
SRC1 = isolhor.C
OBJ1 = $(SRC1:.C=.o)

EXE2 = ih_berlin
SRC2 = ih_berlin.C
OBJ2 = $(SRC2:.C=.o)

$(EXE1): $(OBJ1) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)
$(EXE2): $(OBJ2) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)




.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1)  $(OBJ2) 
	rm -f $(EXE1) $(EXE2)

