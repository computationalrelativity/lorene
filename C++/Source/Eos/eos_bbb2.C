/*
 *  Methods of class Eos_BBB2
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


char eos_bbb2_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2001/09/11  16:23:00  eric
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
Eos_BBB2::Eos_BBB2(const char* path)
		: Eos_tabul(
		   "EOS BBB2 [Baldo, Bombaci & Burgio (1997)]",
		            "eos_bbb2.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_BBB2::Eos_BBB2(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_BBB2::Eos_BBB2(ifstream& fich) : Eos_tabul(fich, "eos_bbb2.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_BBB2::~Eos_BBB2(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_BBB2::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_BBB2 !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_BBB2::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_BBB2::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_BBB2 (Baldo, Bombaci & Burgio 1997) : "
    	<< endl ;
    	
    ost << "  model : BHF (Paris +TBF) " << endl ;

    return ost ;

}

			
