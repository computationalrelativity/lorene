/*
 * Methods for Base_vect and file manipulation
 *
 * (see file base_vect.h for documentation)
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char base_vect_from_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/02/09  13:25:23  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>

// Headers C
#include <stdlib.h>

// Header Lorene
#include "base_vect.h"
#include "utilitaires.h"

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//


int Base_vect_cart::identify() const	{ return 1; }

int Base_vect_spher::identify() const	{ return 2; }



		//--------------------------------------//
		//  Base_vect construction from a file	//
		//--------------------------------------//

Base_vect* Base_vect::bvect_from_file(FILE* fich) {
    
    Base_vect* p_bvect ; 
    
    // Type (class) of vectorial basis identificator ; 
    int identificator ;     
    fread_be(&identificator, sizeof(int), 1, fich) ;		

    switch(identificator) {
	
	case 1 : {
	    p_bvect = new Base_vect_cart(fich) ; 
	    break ; 
	}
	
	case 2 : {
	    p_bvect = new Base_vect_spher(fich) ; 
	    break ; 
	}
	
	default : {
	    cout << "Base_vect::bvect_from_file : unknown type of Base_vect!" 
		 << endl ; 
	    cout << " identificator = " << identificator << endl ; 
	    abort() ; 
	    break ; 
	}
	
    }
    
    return p_bvect ; 
    
}

