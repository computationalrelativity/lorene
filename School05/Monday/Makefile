SHELL=/bin/sh

.SUFFIXES : .o .C

SRC = grid.C grid_uniform.C grid_chebyshev_gauss.C grid_chebyshev_gl.C \
      grid_legendre_gauss.C grid_legendre_gl.C ortho_poly.C chebyshev_poly.C \
      legendre_poly.C coef_cheb_fft.C legendre_nodes.C  plot.C
 
OBJ = $(SRC:.C=.o)

SRC1 = es.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = es

SRC2 = runge.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = runge

SRC3 = cheby.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = cheby

SRC4 = legendre.C
OBJ4 = $(SRC4:.C=.o)
EXE4 = legendre

SRC5 = interpole.C
OBJ5 = $(SRC5:.C=.o)
EXE5 = interpole

SRC6 = approx_ortho.C
OBJ6 = $(SRC6:.C=.o)
EXE6 = approx_ortho

# C++ compiler:
# ------------
CXX 	 = g++

# Options for the C++ compiler for debugging:
# -------------------------------------------
CXXFLAGS_G = -g -pedantic -Wall -W -Wundef -Wshadow -Wcast-qual \
  -Wcast-align -Wconversion -Winline \
  -Wabi -Wold-style-cast -Woverloaded-virtual 

# Path for the include files (e.g. for PGPLOT):
# --------------------------------------------
INC	= -I/usr/local/include

# Graphical libraries: PGPLOT and X11
# -----------------------------------
#
# NB1: change -lg2c to -lgfortran if you are using gcc >= 4.0
# NB2: add -lpng if your version of PGPLOT requires PNG support
#
LIB_GRAPH = -L/usr/X11R6/lib -lcpgplot -lpgplot -lgfortran -lX11

# C, C++ libarary, mathematical library
# -------------------------------------
LIB_CXX = -lstdc++ -lfftw3 -lm


$(EXE1): $(OBJ1) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)

$(EXE2): $(OBJ2) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)

$(EXE3): $(OBJ3) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)

$(EXE4): $(OBJ4) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)

$(EXE5): $(OBJ5) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)

$(EXE6): $(OBJ6) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ6) $(OBJ) $(LIB_GRAPH) $(LIB_CXX)


all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

doc: doc.C $(SRC)  *.h
	doxygen doxyfile

clean:
	rm -f $(OBJ) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) 
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) 
	rm -fr Doc
	rm -fr *~

