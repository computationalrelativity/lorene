/*
 *  Method of class Eos_fit_SLy4
 *
 *    (see file eos_fitting.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char eos_fit_sly4_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/09/26 18:54:35  k_taniguchi
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "headcpp.h"
#include "eos.h"
#include "eos_fitting.h"

//--------------------------------//
//          Constructors          //
//--------------------------------//

// Standard constructor
// --------------------
Eos_fit_SLy4::Eos_fit_SLy4(const char* path)
    : Eos_fitting("EOS fitted to SLy4", "eos_fit_sly4.d", path)
{}

// Constructor from binary file
// ----------------------------
Eos_fit_SLy4::Eos_fit_SLy4(FILE* fich) : Eos_fitting(fich) {}

// Constructor from a formatted file
// ---------------------------------
Eos_fit_SLy4::Eos_fit_SLy4(ifstream& fich)
    : Eos_fitting(fich, "eos_fit_sly4.d")
{}

          //------------------------------//
          //          Destructor          //
          //------------------------------//

Eos_fit_SLy4::~Eos_fit_SLy4() {

    // does nothing

}

          //----------------------------------------//
          //          Comparison operators          //
          //----------------------------------------//

bool Eos_fit_SLy4::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
        cout << "The second EOS is not of type Eos_fit_SLy4 !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_fit_SLy4::operator!=(const Eos& eos_i) const {

  return !(operator==(eos_i)) ;

}

          //---------------------------//
          //          Outputs          //
          //---------------------------//

ostream& Eos_fit_SLy4::operator>>(ostream& ost) const {

    ost <<
      "EOS of class Eos_fit_SLy4 : "
	<< endl ;

    ost << "  composition : n, p, e, mu" << endl ;
    ost << "  model : effective nucleon energy functional, SLy4" << endl ;

    return ost ;

}
