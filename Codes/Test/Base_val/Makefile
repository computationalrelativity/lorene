SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

SRC1 = test_name.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)


$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1) 
	rm -f $(EXE1) 

