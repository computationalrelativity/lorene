/*
 *  Methods for the class Et_bin_bhns_extr
 *
 *    (see file et_bin_bhns_extr.h for documentation).
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

char et_bin_bhns_extr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:48:19  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "et_bin_bhns_extr.h"
#include "etoile.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------
Et_bin_bhns_extr::Et_bin_bhns_extr(Map& mpi, int nzet_i, bool relat,
				   const Eos& eos_i, bool irrot,
				   const Base_vect& ref_triad_i)
  : Etoile_bin(mpi, nzet_i, relat, eos_i, irrot, ref_triad_i)
{}

// Copy constructor
// ----------------
Et_bin_bhns_extr::Et_bin_bhns_extr(const Et_bin_bhns_extr& ns)
  : Etoile_bin(ns)
{}

// Constructor from a file
// -----------------------
Et_bin_bhns_extr::Et_bin_bhns_extr(Map& mpi, const Eos& eos_i,
				   const Base_vect& ref_triad_i, FILE* fich)
  : Etoile_bin(mpi, eos_i, ref_triad_i, fich)
{}

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

}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_bin_bhns_extr::sauve(FILE* fich) const {

    Etoile_bin::sauve(fich) ;

}
