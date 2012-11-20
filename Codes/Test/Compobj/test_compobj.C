/*
 *  Test code for the Compobj and Compobj_QI classes
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
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

char test_compobj_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2012/11/20 16:30:21  c_some
 * Added test for class Compobj_QI
 *
 * Revision 1.1  2012/11/15 16:20:52  c_some
 * New class Compobj
 *
 *
 * $Header$
 *
 */

// C++ headers

// C headers

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"


int main() {
	
	// Number of domains
	int nz = 3 ;

	// Number of coefficients for the spectral expansions (the same in each domain)
	int nr = 9 ; 
	int nt = 7 ; 
	int np = 1 ; 

	// Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // symmetry with respect to phi --> phi + pi
   	bool compact = true ; // external domain is compactified

    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

    cout << mgrid << endl ; 
    
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------
  
  	double r_limits[] = {0, 1., 2., __infinity} ; 
    Map_af map(mgrid, r_limits) ;

    cout << map << endl ; 
    
	// Construction of the compact object : 
	// ----------------------------------
	
	Compobj star(map) ; 
	
	cout << "star :" << star << endl ; 	
	
	Compobj_QI starqi(map) ; 
	
	cout << "starqi :" << starqi << endl ; 	
	
	cout << "ADM mass : " << starqi.mass_g() << endl ;
	
}
