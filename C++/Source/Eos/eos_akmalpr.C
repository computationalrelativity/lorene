/*
 *  Methods of class Eos_AkmalPR
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


char eos_akmalpr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2001/09/11  16:22:41  eric
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
Eos_AkmalPR::Eos_AkmalPR(const char* path)
		: Eos_tabul(
		"EOS AkmalPR [Akmal, Pandharipande & Ravenhall (1998)]",
		            "eos_akmalpr.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_AkmalPR::Eos_AkmalPR(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_AkmalPR::Eos_AkmalPR(ifstream& fich) : Eos_tabul(fich, "eos_akmalpr.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_AkmalPR::~Eos_AkmalPR(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_AkmalPR::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_AkmalPR !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_AkmalPR::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_AkmalPR::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_AkmalPR (Akmal, Pandharipande & Ravenhall 1998) : "
    	<< endl ;
    	
    ost << "  composition :  n,p,e,mu " << endl ;
    ost << "  model : A18+dv+UIX* " << endl ;

    return ost ;

}

			
