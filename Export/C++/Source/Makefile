SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

BIN_TOOLS = $(HOME_LORENE)/bin

SRC	= *.C
INCEXP	  = $(INC) -I$(HOME_LORENE)/Export/C++/Include

install: $(SRC) Objects Objects_g
	touch Dependencies
	makedepend -f Dependencies $(INCEXP) $(SRC)  # determines the dependencies
	echo *C | $(BIN_TOOLS)/pointCpointo > List_obj     # list of *.o files
	mv Objects_g/*.o .
	make -f Makefile_lib_g
	mv *.o Objects_g
	mv Objects/*.o .
	make -f Makefile_lib
	touch $(HOME_LORENE)/Lib/date_mark
	mv *.o Objects

uninstall:
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
