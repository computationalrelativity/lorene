/*
 *  Methods of class Eos_SLy4
 *
 *  (see file eos_tabul.h for documentation).
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


char eos_sly4_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/11/22  19:31:19  eric
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
Eos_SLy4::Eos_SLy4(const char* path)
		: Eos_tabul("EOS SLy4 [Douchin & Haensel (2000)]",
		            "eos_sly4.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_SLy4::Eos_SLy4(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_SLy4::Eos_SLy4(ifstream& fich) : Eos_tabul(fich, "eos_sly4.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_SLy4::~Eos_SLy4(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_SLy4::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_SLy4 !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_SLy4::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_SLy4::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_SLy4 (SLy4 model of Douchin & Haensel 2000) : "
    	<< endl ;
    	
    ost << "  composition :  n,p,e,mu" << endl ;
    ost << "  model : effective nucleon energy functional, SLy4" << endl ;

    return ost ;

}

			
