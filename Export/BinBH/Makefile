SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

INCEXP	  = -I$(HOME_LORENE)/Export/C++/Include
LIB_G = -L$(HOME_LORENE)/Lib -llorene_export_g -llorene_g -llorenef77_g


.SUFFIXES : .o .C

EXE = readinit visu_1d

SRC1 = readinit.C
OBJ1 = $(SRC1:.C=.o)

SRC2 = visu_1d.C
OBJ2 = $(SRC2:.C=.o)

readinit: $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
        $(LIB_GSL) $(LIB_CXX)

visu_1d: $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
        $(LIB_GSL) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INCEXP) $<

uninstall:
	rm -f $(OBJ1) $(OBJ2)
	rm -f $(EXE)

