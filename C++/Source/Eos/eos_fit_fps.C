/*
 *  Method of class Eos_fit_FPS
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:52:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2004/10/24 19:14:53  k_taniguchi
 * Correction of the file name which is called in the constructor from file.
 *
 * Revision 1.1  2004/09/26 18:55:10  k_taniguchi
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
namespace Lorene {
Eos_fit_FPS::Eos_fit_FPS(const char* path)
    : Eos_fitting("EOS fitted to FPS", "eos_fit_fps.d", path)
{}

// Constructor from binary file
// ----------------------------
Eos_fit_FPS::Eos_fit_FPS(FILE* fich) : Eos_fitting(fich) {}

// Constructor from a formatted file
// ---------------------------------
Eos_fit_FPS::Eos_fit_FPS(ifstream& fich)
    : Eos_fitting(fich, "eos_fit_fps.d")
{}

          //------------------------------//
          //          Destructor          //
          //------------------------------//

Eos_fit_FPS::~Eos_fit_FPS() {

    // does nothing

}

          //----------------------------------------//
          //          Comparison operators          //
          //----------------------------------------//

bool Eos_fit_FPS::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
        cout << "The second EOS is not of type Eos_fit_FPS !" << endl ;
	resu = false ;
    }

    return resu ;

}

bool Eos_fit_FPS::operator!=(const Eos& eos_i) const {

  return !(operator==(eos_i)) ;

}

          //---------------------------//
          //          Outputs          //
          //---------------------------//

ostream& Eos_fit_FPS::operator>>(ostream& ost) const {

    ost <<
      "EOS of class Eos_fit_FPS : "
	<< endl ;

    ost << "  composition : n, p, e, mu" << endl ;
    ost << "  model : effective nucleon energy functional, FPS" << endl ;

    return ost ;

}
}
