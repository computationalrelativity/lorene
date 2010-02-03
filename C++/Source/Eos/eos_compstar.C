/*
 *  Methods of class Eos_Compstar
 *
 *  (see file eos_tabul.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010 Jerome Novak
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


char eos_compstar_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/02/03 14:56:45  j_novak
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "headcpp.h"
#include "eos.h"

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
Eos_Compstar::Eos_Compstar(const char* file_name)
		: Eos_tabul("EOS Compstar", file_name)
{}


// Constructor from binary file
// ----------------------------
Eos_Compstar::Eos_Compstar(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_Compstar::Eos_Compstar(ifstream& fich) : Eos_tabul(fich)
 {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_Compstar::~Eos_Compstar(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_Compstar::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_Compstar !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_Compstar::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_Compstar::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_Compstar : tabulated EOS for the 2010 CompStar school."
    	<< endl ;
    	
    return ost ;

}

			
