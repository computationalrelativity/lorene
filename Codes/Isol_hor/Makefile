SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a

.SUFFIXES : .o .C

EXE = isolhor
SRC = isolhor.C
OBJ = $(SRC:.C=.o)

$(EXE): $(OBJ) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)  $(LIB_GSL)




.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ) 
	rm -f $(EXE)

