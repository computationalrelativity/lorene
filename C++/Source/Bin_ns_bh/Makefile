SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

BIN_TOOLS = $(HOME_LORENE)/bin

SRC	= *.C

install: $(SRC) Objects Objects_g
	touch Dependencies
	makedepend -f Dependencies $(INC) $(SRC)  # determines the dependencies
	echo *C | $(BIN_TOOLS)/pointCpointo > List_obj     # list of *.o files
	mv Objects_g/*.o .
	mv _rien.o Objects_g
	make -f Makefile_lib_g
	mv *.o Objects_g
	mv Objects/*.o .
	mv _rien.o Objects
	make -f Makefile_lib
	mv *.o Objects

uninstall:
	rm -f *.o
	rm Dependencies
	touch Dependencies.bak
	rm Dependencies.bak
	rm List_obj
	rm -fr Objects	
	rm -fr Objects_g

Objects:
	mkdir Objects
	touch Objects/_rien.o

Objects_g:
	mkdir Objects_g
	touch Objects_g/_rien.o
