/*
 * Method of class Scalar to check if a Poisson equation has been correctly
 * solved.
 *
 * (see file cmp.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 2000-2001 Eric Gourgoulhon (for a preceding Cmp version)
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


char scalar_test_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/25 09:13:11  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tensor.h"

Tbl Scalar::test_poisson(const Scalar& uu, ostream& ostr, bool detail) const {
    
    assert( &(uu.get_mp()) == mp ) ; 
    
    // Computation of the absolute and relative errors
    // -----------------------------------------------
    
    int dzi ;
    if ( check_dzpuis(4) ) {
	dzi = 4 ; 
    }
    else{
	if ( check_dzpuis(2) ) {
	    dzi = 2 ; 
	}
	else{
	    assert( check_dzpuis(0) ) ; 
	    dzi = 0 ; 
	}
    }
    
    Tbl tdiff = max( abs( uu.laplacien(dzi) - *this ) ) ;

    Tbl tmax = max( abs(*this) ) ; 
    
    int nz = mp->get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ; 
    
    Tbl trel(nz) ; 
    trel.set_etat_qcq() ; 

    if ( (dzpuis == 0) || (tmax(nzm1) == double(0)) ) {

	double s_max = max( tmax ) ; 

	for (int l=0; l<nz; l++) {
	    trel.set(l) = tdiff(l) / s_max ;     
	}

    }
    else{

	double s_max = 0 ; 
	for (int l=0; l<nzm1; l++) {
	    s_max = (tmax(l) > s_max) ? tmax(l) : s_max ;     
	}

	for (int l=0; l<nzm1; l++) {
	    trel.set(l) = tdiff(l) / s_max ;     
	}

	trel.set(nzm1) = tdiff(nzm1) / tmax(nzm1) ; 

    }
    
    // All the errors are set in the same output 2-D Tbl
    // -------------------------------------------------

    Tbl err(3, nz) ; 
    
    err.set_etat_qcq() ; 
    
    for(int l=0; l<nz; l++) {
	err.set(0, l) = trel(l) ; 
	err.set(1, l) = tdiff(l) ; 
	err.set(2, l) = tmax(l) ; 
    }

    // Display
    // -------
    
    if (detail) {
	ostr << "Max. source :" ;
	for (int l=0; l<nz; l++) {
	    ostr << "  " << err(2, l) ;
	}
    
	ostr << endl << "Abs. error : " ;
	for (int l=0; l<nz; l++) {
	    ostr << "  " << err(1, l) ;
	}
    }
    
    ostr << endl << "Rel. error : " ;
    for (int l=0; l<nz; l++) {
	ostr << "  " << err(0, l) ;
    }
    
    ostr << endl ; 
    
    return err ; 

}

