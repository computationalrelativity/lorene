SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

INCEXP	  = -I$(HOME_LORENE)/Export/C++/Include
LIB_G = -L$(HOME_LORENE)/Lib -llorene_export_g -llorene_g -llorenef77_g


.SUFFIXES : .o .C

EXE = read_bin_ns
SRC = read_bin_ns.C
OBJ = $(SRC:.C=.o)

$(EXE): $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INCEXP) $<

uninstall:
	rm -f $(OBJ)
	rm -f $(EXE)

