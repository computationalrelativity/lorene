/*
 *  Methods of class Base_vect
 *
 *   (see file bse_vect.h for documentation)
 *
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


char base_vect_C[] = "$Header$" ;

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
 * Revision 2.2  2000/02/09  13:24:12  eric
 * REFONTE COMPLETE DE LA CLASSE
 * L'identification n'est plus base sur un membre statique (numero
 * d'instance) mais sur les caracteres physiques (rot_phi, etc...)
 * Ajout des constructeurs par copie et lecture de fichier.
 *
 * Revision 2.1  2000/01/10  15:43:11  eric
 * Methode change_basis (bidon).
 *
 * Revision 2.0  2000/01/10  12:43:21  eric
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
#include <string.h>

// Headers Lorene
#include "base_vect.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor without name
// ---------------------------------
Base_vect::Base_vect(){
    
    set_name("") ; 
    
}

// Standard constructor with name
// ------------------------------
Base_vect::Base_vect(const char* name_i){
    
    set_name(name_i) ; 
    
}


// Copy constructor
// ----------------
Base_vect::Base_vect(const Base_vect& bvect_i){
        
    set_name(bvect_i.name) ; 
    
}

// Constructor from file
// ---------------------
Base_vect::Base_vect(FILE* fich){
        
    fread(name, sizeof(char), 100, fich) ;		
    
}


			//--------------//
			//  Destructor  //
			//--------------//

Base_vect::~Base_vect(){
    
    // does nothing
        
}

			//-------------------------//
			//  Manipulation of name   //
			//-------------------------//
			
			
void Base_vect::set_name(const char* name_i) {

    strncpy(name, name_i,  100) ; 
    
}

const char* Base_vect::get_name() const {
    
    return name ; 
    
}

			//------------//
			//  Outputs   //
			//------------//

void Base_vect::sauve(FILE* fich) const {

    int ident = identify() ; 
    fwrite_be(&ident, sizeof(int), 1, fich) ;	
    	
    fwrite(name, sizeof(char), 100, fich) ;		
   
}
    



ostream& operator<<(ostream& ost, const Base_vect& bvect)  {
    ost << bvect.get_name() << endl ; 
    bvect >> ost ;
    return ost ;
}



		    //----------------------//
		    // Comparison operator  //
		    //----------------------//
		    
bool Base_vect::operator!=(const Base_vect& bi) const {
    
    return !(bi == *this) ; 
    
}

