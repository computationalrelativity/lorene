/*
 *  Code for testing the covariant derivatives through the Connection class.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char test_connect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/03 14:27:51  e_gourgoulhon
 * First non trivial test (successfull !).
 *
 * Revision 1.1  2003/10/02 21:33:02  e_gourgoulhon
 * Test code for Connection.
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "connection.h"
#include "nbr_spx.h"
#include "utilitaires.h"


int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nr = 5 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 4 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1., 2., __infinity} ; 
  	assert( nz == 3 ) ;  // since the above array described only 3 domains
  
	Map_af map(mgrid, r_limits) ; 
  
  
	// Construction of a flat connection
	// ---------------------------------
	
	// Representation on a spherical orthonormal basis
	Connection_fspher ders(map, map.get_bvect_spher()) ; 
	
	// Representation on a Cartesian orthonormal basis
	Connection_fcart derc(map, map.get_bvect_cart()) ; 


	// Construction of a scalar field (Scalar)
	// ---------------------------------------

	const Coord& x = map.x ; 

	Scalar uu(map) ; 

	uu = x ; 	
	uu.annule(nz-1,nz-1) ; 	// zero in the last domain
	
	uu.std_spectral_base_scal() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	cout << "uu : " << uu << endl ; 
	uu.spectral_display(cout) ; 
	arrete() ; 
	
	// Gradient of the scalar field
	// ----------------------------
	
	Vector duuc = derc.derive_cov(uu)  ;  
	cout << "duuc : " << duuc << endl ; 

	arrete() ; 

	Vector duus = ders.derive_cov(uu) ; 
	cout << "duus : " << duus << endl ; 

	arrete() ; 
	
	// Test
	// ----
	
	Vector duus_c = duus ; 
	duus_c.change_triad( map.get_bvect_cart() ) ; 
	
	Vector diffc = duus_c - duuc ; 
	
	cout << "diffc : " << diffc << endl ; 
	cout << "Norm of diffc: " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme(diffc(i)) << endl ; 
	}
	arrete() ; 
	
	
	Vector duuc_s = duuc ; 
	duuc_s.change_triad( map.get_bvect_spher() ) ; 
	
	Vector diffs = duuc_s - duus ; 

	cout << "diffs : " << diffs << endl ; 
	cout << "Norm of diffs: " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme(diffs(i)) << endl ; 
	}
	

	return EXIT_SUCCESS ; 
}
