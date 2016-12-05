/*
 *  Methods for the class Et_bin_bhns_extr
 *
 *    (see file et_bin_bhns_extr.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004-2005 Keisuke Taniguchi
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
 * Revision 1.5  2016/12/05 16:17:52  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:54  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:07  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/02/28 23:09:38  k_taniguchi
 * Modification of some functions to include the case of the conformally flat
 * background metric.
 *
 * Revision 1.1  2004/11/30 20:48:19  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Et_bin_bhns_extr::Et_bin_bhns_extr(Map& mpi, int nzet_i, bool relat,
				   const Eos& eos_i, bool irrot,
				   const Base_vect& ref_triad_i,
				   bool kerrs, bool multi)
  : Etoile_bin(mpi, nzet_i, relat, eos_i, irrot, ref_triad_i),
    kerrschild(kerrs),
    multipole(multi)
{}

// Copy constructor
// ----------------
Et_bin_bhns_extr::Et_bin_bhns_extr(const Et_bin_bhns_extr& ns)
  : Etoile_bin(ns),
    kerrschild(ns.kerrschild),
    multipole(ns.multipole)
{}

// Constructor from a file
// -----------------------
Et_bin_bhns_extr::Et_bin_bhns_extr(Map& mpi, const Eos& eos_i,
				   const Base_vect& ref_triad_i, FILE* fich)
  : Etoile_bin(mpi, eos_i, ref_triad_i, fich)
{

    // kerrschild is read in the file:
    fread(&kerrschild, sizeof(bool), 1, fich) ;

    // multipole is read in the file:
    fread(&multipole, sizeof(bool), 1, fich) ;

}

			    //------------//
			    // Destructor //
			    //------------//

Et_bin_bhns_extr::~Et_bin_bhns_extr()
{}

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_bin_bhns_extr
// --------------------------------------
void Et_bin_bhns_extr::operator=(const Et_bin_bhns_extr& ns) {

    // Assignment of quantities common to the derived classes of Etoile_bin
    Etoile_bin::operator=(ns) ;

    // Assignment of proper quantities of class Et_bin_bhns_extr
    kerrschild = ns.kerrschild ;
    multipole = ns.multipole ;

}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_bin_bhns_extr::sauve(FILE* fich) const {

    Etoile_bin::sauve(fich) ;

    fwrite(&kerrschild, sizeof(bool), 1, fich) ;

    fwrite(&multipole, sizeof(bool), 1, fich) ;

}
}
