/*
 *  Methods of class Eos_BalbN1H1
 *
 *  (see file eos_tabul.h for documentation).
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
 *   Copyright (c) 2001 J. Leszek Zdunik
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


char eos_balbn1h1_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2001/09/11  16:22:50  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C++
#include <iostream.h>

// Headers Lorene
#include "eos.h"

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_BalbN1H1::Eos_BalbN1H1(const char* path)
		: Eos_tabul(
		"EOS BalbN1H1 [Balberg  (2000)]",
		            "eos_balbn1h1.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_BalbN1H1::Eos_BalbN1H1(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_BalbN1H1::Eos_BalbN1H1(ifstream& fich) : Eos_tabul(fich, "eos_balbn1h1.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_BalbN1H1::~Eos_BalbN1H1(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_BalbN1H1::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_BalbN1H1 !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_BalbN1H1::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_BalbN1H1::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_BalbN1H1 (Balberg 2000) : "
    	<< endl ;
    	
    return ost ;

}

			
