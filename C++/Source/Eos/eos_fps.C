/*
 *  Methods of class Eos_FPS
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


char eos_fps_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/11/22  19:30:29  eric
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
namespace Lorene {
Eos_FPS::Eos_FPS(const char* path)
		: Eos_tabul("EOS FPS",
		            "eos_fps.d", path)
{}


// Constructor from binary file
// ----------------------------
Eos_FPS::Eos_FPS(FILE* fich) : Eos_tabul(fich) {}



// Constructor from a formatted file
// ---------------------------------
Eos_FPS::Eos_FPS(ifstream& fich) : Eos_tabul(fich, "eos_fps.d") {}



			//--------------//
			//  Destructor  //
			//--------------//

Eos_FPS::~Eos_FPS(){

    // does nothing

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_FPS::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_FPS !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_FPS::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}

			//------------//
			//  Outputs   //
			//------------//


ostream& Eos_FPS::operator>>(ostream & ost) const {

    ost <<
    "EOS of class Eos_FPS : Friedman-Pandharipande + Skyrme "
    	<< endl ;
    ost << "  (Lorenz, Ravenhall, Petick)" << endl ;
    	
    ost << "  composition :  n,p,e,mu" << endl ;
    ost << "  model : effective nucleon energy functional, FPS" << endl ;
    ost << "  BPS EOS below neutron drip point" << endl ;
    ost << "  Crust bottom at n = 0.0957 fm^{-3}" << endl ;

    return ost ;

}

			
}
